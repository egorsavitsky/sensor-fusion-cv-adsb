#pragma once

#include "fusion.hpp"
#include <ostream>
#include <string>

std::string trackCsvHeader(const std::string& prefix);
void writeTrackCsv(std::ostream& out, const Track3D* track);
const Track3D& findTrack(const std::vector<Track3D>& tracks, int id);
void printOverlapPair(std::ostream& out, double time, const Track3D& cv,
                      const Track3D& adsb, const OverlapPair& pair);
void writeOverlapPair(std::ostream& out, double time, const Track3D& cv,
                      const Track3D& adsb, const OverlapPair& pair, bool correct);
void writeAssociationRow(std::ostream& out, double time, const char* status,
                         const Track3D* cv, const Track3D* adsb, const Match* match, bool correct);
