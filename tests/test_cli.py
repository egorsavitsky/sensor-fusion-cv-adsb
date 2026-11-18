"""Проверка терминального отчёта и CSV на примере с самолётом и машиной.

Запуск: python tests/test_cli.py build/fusion_core.exe
"""

import csv
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

repo = Path(__file__).resolve().parent.parent
program = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else repo / 'build' / ('fusion_core.exe' if os.name == 'nt' else 'fusion_core')
assert program.is_file(), f'Build the program first: {program}'
build = (repo / 'build').resolve()
assert build.parent == repo, 'Temporary test output must stay in the repository'
# Шум АЗН-В отключён, чтобы точные числа пар не зависели от реализации normal_distribution.
base = [str(program), '--demo', '--assignment', 'optimal', '--horizontal-sigma', '0',
        '--vertical-sigma', '0', '--heading-sigma', '0', '--dropout', '0', '--coverage', '1']

# Все выходные файлы создаются в отдельной папке внутри build и удаляются после проверки.
with tempfile.TemporaryDirectory(prefix='overlap-cli-', dir=build) as directory:
    temporary = Path(directory)
    results = []
    for index, threshold in enumerate([0, .1, .8]):
        paths = {name: temporary / f'{index}-{name}.csv' for name in ['summary', 'pairs', 'matches', 'model-log', 'cv-log']}
        arguments = base + ['--overlap-threshold', str(threshold)]
        for name, path in paths.items():
            arguments += ['--' + name, str(path)]
        process = subprocess.run(arguments, cwd=repo, capture_output=True, text=True, encoding='utf-8', check=True)
        tables = {}
        for name, path in paths.items():
            with path.open(encoding='utf-8', newline='') as stream:
                tables[name] = list(csv.DictReader(stream))
            assert all(None not in row and all(value is not None for value in row.values()) for row in tables[name]), name
        assert len(tables['summary']) == 1
        stats = tables['summary'][0]
        rows = tables['pairs']
        matched = [row for row in tables['matches'] if row['status'] == 'matched']
        results.append(stats)

        # Счётчики в итогах сверяются с полными таблицами, а не только с текстом консоли.
        assert int(stats['evaluation_frames']) == 12 and int(stats['observations']) == 24
        assert int(stats['synthetic_objects']) == 2 and int(stats['unique_cv']) == 2 and int(stats['unique_adsb']) == 2
        assert int(stats['positive_pairs']) == len(rows)
        assert int(stats['threshold_pairs']) == sum(row['above_threshold'] == '1' for row in rows)
        assert int(stats['eligible_pairs']) == sum(row['above_threshold'] == row['same_class'] == row['height_overlap'] == '1' for row in rows)
        assert int(stats['matches']) == len(matched) == sum(row['assigned'] == '1' for row in rows)
        assert int(stats['correct']) + int(stats['wrong']) == len(matched)
        assert int(stats['matches']) + int(stats['cv_only']) == int(stats['observations'])
        assert int(stats['matches']) + int(stats['adsb_only']) == int(stats['adsb_observations'])
        assert stats['camera_right_m'] == '0' and stats['camera_up_m'] == '1.1' and stats['camera_forward_m'] == '15.5'
        assert float(stats['calibration_frames']) == 0
        assert len(tables['cv-log']) == 24 and len(tables['model-log']) == int(stats['synthetic_messages'])
        assert {row['class'] for row in tables['model-log']} == {'0', '1'}
        assert len({(row['timestamp_s'], row['cv_id']) for row in matched}) == len(matched)
        assert len({(row['timestamp_s'], row['adsb_id']) for row in matched}) == len(matched)
        for row in rows:
            assert 0 < float(row['iou']) <= 1 and float(row['intersection_m2']) > 0
            assert row['overlap_positive'] == '1'
            assert (row['above_threshold'] == '1') == (float(row['iou']) >= threshold)
            for prefix in ['cv_', 'adsb_']:
                assert all(float(row[prefix + dimension + '_m']) > 0 for dimension in ['width', 'height', 'length'])
                assert math.isfinite(float(row[prefix + 'heading_deg']))
        for row in tables['model-log']:
            assert float(row['received_at_s']) >= float(row['measurement_time_s'])

        printed = [line for line in process.stdout.splitlines() if line.startswith('t=')]
        assert len(printed) == int(stats['threshold_pairs'])
        for line, row in zip(printed, (r for r in rows if r['above_threshold'] == '1')):
            time, cv, adsb, scores = re.fullmatch(r't=(\S+) CV=\{([^}]+)\} ADS-B=\{([^}]+)\} (.*)', line).groups()
            assert float(time) == float(row['timestamp_s'])
            for prefix, track in [('cv_', cv), ('adsb_', adsb)]:
                identifier, cls, label, pos, size, heading = re.fullmatch(
                    r'id=(-?\d+), class=(-?\d+)\((\w+)\), pos=\(([^)]+)\), size=\(([^)]+)\), heading_deg=(\S+)', track).groups()
                assert (identifier, cls) == (row[prefix + 'id'], row[prefix + 'class'])
                assert label == {'0': 'aircraft', '1': 'vehicle'}[cls]
                values = pos.split(',') + size.split(',') + [heading]
                fields = ['east_m', 'north_m', 'up_m', 'width_m', 'height_m', 'length_m', 'heading_deg']
                assert len(values) == len(fields)
                for value, field in zip(values, fields):
                    assert math.isclose(float(value), float(row[prefix + field]), rel_tol=6e-6, abs_tol=1e-12)
            flags = dict(item.split('=') for item in scores.split())
            for field in ['overlap_positive', 'above_threshold', 'height_overlap', 'same_class', 'assigned']:
                assert flags[field] == row[field]
            for field in ['iou', 'intersection_m2']:
                assert math.isclose(float(flags[field]), float(row[field]), rel_tol=6e-6)
        assert f"Evaluated frames: {stats['evaluation_frames']}, CV observations: {stats['observations']}" in process.stdout
        assert (f"CV/ADS-B observations: {stats['observations']} / {stats['adsb_observations']}, "
                f"unique IDs: {stats['unique_cv']} / {stats['unique_adsb']}") in process.stdout
        assert (f"Positive overlap pairs: {len(rows)}, at threshold: {stats['threshold_pairs']}, "
                f"eligible: {stats['eligible_pairs']}") in process.stdout
        assert f"Matches: {len(matched)}, correct: {stats['correct']}, wrong: {stats['wrong']}" in process.stdout

    assert int(results[0]['threshold_pairs']) >= int(results[1]['threshold_pairs']) >= int(results[2]['threshold_pairs'])
    assert len({row['positive_pairs'] for row in results}) == 1
    assert len({row['synthetic_messages'] for row in results}) == 1
    assert int(results[1]['correct']) == int(results[1]['matches']) == 22

    # Проверяются альтернативный вывод, офлайн-режим и отсутствие второго источника.
    all_pairs = subprocess.run(base + ['--overlap-threshold', '.8', '--all-overlaps'], cwd=repo,
                               capture_output=True, text=True, encoding='utf-8', check=True)
    assert sum(line.startswith('t=') for line in all_pairs.stdout.splitlines()) == int(results[2]['positive_pairs'])
    for name, extra in [('offline', ['--offline']), ('no-coverage', ['--coverage', '0']), ('loss', ['--dropout', '1'])]:
        path = temporary / (name + '.csv')
        process = subprocess.run(base + extra + ['--quiet', '--summary', str(path)], cwd=repo,
                                 capture_output=True, text=True, encoding='utf-8', check=True)
        with path.open() as stream:
            row = next(csv.DictReader(stream))
        assert not any(line.startswith('t=') for line in process.stdout.splitlines())
        if name == 'offline':
            assert int(row['adsb_observations']) == 24
        else:
            assert int(row['matches']) == int(row['positive_pairs']) == int(row['adsb_observations']) == 0
            assert int(row['cv_only']) == 24 and math.isnan(float(row['precision']))

    # Ошибки параметров и совпадение путей не должны приводить к перезаписи исходных файлов.
    fixture = repo / 'tests/data/vertical.tsv'
    original = fixture.read_bytes()
    for extra in [['--overlap-threshold', '-0.1'], ['--overlap-threshold', '1.1'], ['--overlap-threshold', 'nan'],
                  ['--period', '0'], ['--input', str(fixture)], ['--pairs', str(fixture)]]:
        failure = subprocess.run(base + extra, cwd=repo, capture_output=True)
        assert failure.returncode != 0, extra
    assert fixture.read_bytes() == original

    alias = temporary / 'shared.csv'
    alias_with_dot = str(temporary) + os.sep + '.' + os.sep + alias.name
    failure = subprocess.run(base + ['--summary', str(alias), '--pairs', alias_with_dot], cwd=repo, capture_output=True)
    assert failure.returncode != 0 and not alias.exists()

    # Пара соприкасается почти вершинами. Положительная площадь должна читаться и из CSV, и из консоли.
    with fixture.open() as stream:
        obstacle = next(csv.DictReader(stream, delimiter='\t'))
    obstacle.update({'Position.X': '1.9999999', 'Position.Z': '1.9999999'})
    corner = temporary / 'corner.tsv'
    with corner.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=obstacle.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(obstacle)
    pairs = temporary / 'corner.csv'
    args = [arg for arg in base if arg != '--demo'] + ['--input', str(corner), '--pairs', str(pairs),
            '--camera-up', '0', '--camera-forward', '0', '--delay', '0', '--overlap-threshold', '0']
    output = subprocess.run(args, cwd=repo, capture_output=True, text=True, check=True).stdout
    with pairs.open() as stream:
        row, = csv.DictReader(stream)
    for field in ['intersection_m2', 'iou']:
        value = float(re.search(r' ' + field + r'=(\S+)', output).group(1))
        assert 0 < value < 1e-12 and math.isclose(value, float(row[field]), rel_tol=6e-6)
    assert row['assigned'] == '1'

    obstacle['ID'] = ''
    with corner.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=obstacle.keys(), delimiter='\t')
        writer.writeheader()
        writer.writerow(obstacle)
    broken = subprocess.run([str(program), '--input', str(corner)], cwd=repo, capture_output=True, text=True)
    assert broken.returncode != 0 and 'Log line 2:' in broken.stderr

print('CLI checks passed: thresholds 0/0.1/0.8, terminal fields, CSV totals, shape metadata, delays, loss and output protection')
