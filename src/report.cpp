#include "report.hpp"

#include <stdexcept>

namespace {
/** Возвращает подпись класса для терминала.
 * @param object_class 0 — самолёт, 1 — машина, остальные значения — неизвестный класс.
 * @return Краткое название для терминала.
 */
const char* className(int object_class) {
    if (object_class == 0) return "aircraft";
    if (object_class == 1) return "vehicle";
    return "unknown";
}

/** Выводит атрибуты одного объекта в читаемом виде.
 * @param out Поток терминального отчёта.
 * @param track Объект в ENU; размер задаётся как ширина, высота, длина.
 */
void printTrack(std::ostream& out, const Track3D& track) {
    out << "{id=" << track.id << ", class=" << track.object_class << '(' << className(track.object_class)
        << "), pos=(" << track.position.x << ',' << track.position.y << ',' << track.position.z
        << "), size=(" << track.size.x << ',' << track.size.y << ',' << track.size.z
        << "), heading_deg=" << track.heading*180/3.14159265358979323846 << '}';
}
} // namespace

/** Формирует названия столбцов одного объекта.
 * @param prefix Префикс источника: cv_, adsb_ либо пустая строка.
 * @return Девять разделённых запятыми названий, без конечной запятой.
 */
std::string trackCsvHeader(const std::string& prefix) {
    std::string result;
    for (const auto* name : {"id","class","east_m","north_m","up_m","width_m","height_m","length_m","heading_deg"}) {
        if (!result.empty()) result += ',';
        result += prefix+name;
    }
    return result;
}

/** Записывает координаты, размеры, курс и класс объекта.
 * @param out Выходной CSV-поток с настроенной точностью.
 * @param track Объект либо nullptr для девяти пустых полей.
 */
void writeTrackCsv(std::ostream& out, const Track3D* track) {
    if (!track) { out << ",,,,,,,,"; return; }
    out << track->id << ',' << track->object_class << ',' << track->position.x << ',' << track->position.y
        << ',' << track->position.z << ',' << track->size.x << ',' << track->size.y << ',' << track->size.z
        << ',' << track->heading*180/3.14159265358979323846;
}

/** Находит объект для записи подробного отчёта.
 * @param tracks Объекты одного источника на одном кадре.
 * @param id Идентификатор искомого объекта.
 * @return Ссылка на объект; при отсутствии ID выбрасывается logic_error.
 */
const Track3D& findTrack(const std::vector<Track3D>& tracks, int id) {
    for (const auto& track : tracks) if (track.id == id) return track;
    throw std::logic_error("Association refers to an unknown track");
}

/** Выводит в терминал обнаруженное пересечение и решение алгоритма.
 * @param out Поток отчёта.
 * @param time Общая временная метка в секундах.
 * @param cv Объект СТЗ.
 * @param adsb Согласованный по времени объект АЗН-В.
 * @param pair Площадь, IoU, признаки фильтрации и назначения.
 */
void printOverlapPair(std::ostream& out, double time, const Track3D& cv,
                      const Track3D& adsb, const OverlapPair& pair) {
    out << "t=" << time << " CV=";
    printTrack(out,cv);
    out << " ADS-B=";
    printTrack(out,adsb);
    out << " intersection_m2=" << pair.intersection_area << " iou=" << pair.iou
        << " overlap_positive=" << (pair.iou > 0) << " above_threshold=" << pair.above_threshold
        << " height_overlap=" << pair.height_overlap << " same_class=" << pair.same_class
        << " assigned=" << pair.assigned << '\n';
}

/** Сохраняет одну пару с положительной площадью пересечения.
 * @param out CSV-поток таблицы всех пересечений.
 * @param time Время кадра.
 * @param cv Объект СТЗ в ENU.
 * @param adsb Объект АЗН-В в ENU.
 * @param pair Геометрия и решение о назначении.
 * @param correct Соответствие опорным ID; используется только в оценке модели.
 */
void writeOverlapPair(std::ostream& out, double time, const Track3D& cv,
                      const Track3D& adsb, const OverlapPair& pair, bool correct) {
    out << time << ',';
    writeTrackCsv(out,&cv);
    out << ',';
    writeTrackCsv(out,&adsb);
    out << ',' << pair.intersection_area << ',' << pair.iou << ',' << (pair.iou > 0)
        << ',' << pair.above_threshold << ',' << pair.height_overlap << ',' << pair.same_class
        << ',' << pair.assigned << ',' << correct << '\n';
}

/** Сохраняет выбранную пару либо объект без пары.
 * @param out CSV-поток результатов назначения.
 * @param time Время кадра.
 * @param status matched, cv_only либо adsb_only.
 * @param cv Объект СТЗ или nullptr.
 * @param adsb Объект АЗН-В или nullptr.
 * @param match Выбранная пара либо nullptr для одиночного объекта.
 * @param correct Правильность назначения по опорным ID модели.
 */
void writeAssociationRow(std::ostream& out, double time, const char* status,
                         const Track3D* cv, const Track3D* adsb, const Match* match, bool correct) {
    out << time << ',' << status << ',';
    writeTrackCsv(out,cv);
    out << ',';
    writeTrackCsv(out,adsb);
    if (match) out << ',' << match->distance << ',' << match->score << ',' << match->overlap
                   << ',' << match->intersection_area << ',' << correct;
    else out << ",,,,,";
    out << '\n';
}
