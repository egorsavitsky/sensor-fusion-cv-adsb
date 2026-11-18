# Данные

`--input` читает TSV с заголовками в формате `URSS_Obstacles_2.txt`. Полная запись СТЗ хранится отдельно от репозитория. [Поля и единицы измерения](../docs/requirements.md#входные-данные-и-координаты)

Эти таблицы получены на примере с самолётом и машиной при seed 42, пороге 0,1 и глобальном назначении.

| Файл | Содержимое |
| --- | --- |
| [demo_cv.csv](demo_cv.csv) | 24 наблюдения СТЗ после перевода в ENU. |
| [demo_adsb.csv](demo_adsb.csv) | 12 модельных сообщений АЗН-В с временем измерения и доставки. |
| [demo_pairs.csv](demo_pairs.csv) | 24 положительных пересечения, из них назначены 22 пары. |

Из папки проекта выгрузки можно получить одной командой.

```powershell
.\build\fusion_core.exe --demo --assignment optimal --quiet --cv-log build\cv.csv --model-log build\adsb.csv --pairs build\pairs.csv
```
