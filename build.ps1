param([switch]$Test, [switch]$Debug)

$ErrorActionPreference = 'Stop'
Push-Location $PSScriptRoot
try {
    New-Item -ItemType Directory -Path build -Force | Out-Null
    $flags = @('-std=c++17', '-Wall', '-Wextra', '-Wpedantic', '-Iinclude')
    if ($Debug) { $flags += @('-O0', '-g') } else { $flags += '-O2' }
    $sources = @('src/fusion.cpp', 'src/geometry.cpp', 'src/log_reader.cpp', 'src/simulation.cpp')
    & g++ @flags @sources src/main.cpp -o build/fusion_core.exe
    if ($LASTEXITCODE -ne 0) { throw 'Application build failed' }
    if ($Test) {
        & g++ @flags @sources tests/test_fusion.cpp -o build/fusion_tests.exe
        if ($LASTEXITCODE -ne 0) { throw 'Test build failed' }
        & ./build/fusion_tests.exe
        if ($LASTEXITCODE -ne 0) { throw 'Tests failed' }
        $result = & ./build/fusion_core.exe --input tests/data/vertical.tsv --calibration-frames 0 --horizontal-sigma 0 --vertical-sigma 0 --coverage 1 --dropout 0 --delay 0
        if ($LASTEXITCODE -ne 0 -or $result -notcontains 'Matches: 6, correct: 6, wrong: 0') {
            throw 'Vertical model integration test failed'
        }
        Write-Output 'Vertical model integration test passed'
    }
} finally {
    Pop-Location
}
