/******************************************************************************
   Copyright (C) 2015-2018 Einar J.M. Baumann <einar.baumann@gmail.com>

   This file is part of the FieldOpt project.

   FieldOpt is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   FieldOpt is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FieldOpt.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#include "welsegs.h"

namespace Simulation {
namespace ECLDriverParts {
Welsegs::Welsegs(Well *well) {
    head_ = "WELSEGS\n";
    foot_ = "/\n\n";
    WelsegsKeyword kw;
    kw.heel_entry = createHeelEntry(well);
    auto segs = well->GetSegments();
    for (int i = 1; i < segs.size(); ++i) {
        kw.seg_entries.push_back(createSegmentEntry(segs[i]));
    }
    keywords_.push_back(kw);
}
QString Welsegs::GetPartString() const {
    if (keywords_.size() == 0)
        return "";

    QString all_keywords = "";
    for (auto kw : keywords_) {
        all_keywords += kw.buildKeyword();
    }
    return all_keywords;
}
QString Welsegs::createHeelEntry(Well *well) {
 /*
  * 0.   Name of well
  * 1.   TVD of top segment.
  * 2. * MD to first segment.
  * 3.   Effective wellbore volume of top segment.
  * 4.   Type of tubing and depth information: INC(remental) or ABS(olute).
  * 5.   Components of pressure drop: 'HFA', 'HF-' or 'H--'.
  * 6.   Flow model: HO or DF.
  * 7.   X-coordinate of nodal point.
  * 8.   Y-coordinate of nodal point.
  */
    auto rseg = well->GetSegments()[0];
    auto record = GetBaseEntryLine(10);
    record[0] = well->name();
    record[1] = QString::number(rseg.TVDChange());
    record[4] = "INC";
    record[5] = "HF-";
    record[6] = "DF";
    record[7] = QString::number(well->trajectory()->GetWellBlocks()->at(0)->getEntryPoint().x());
    record[8] = QString::number(well->trajectory()->GetWellBlocks()->at(0)->getEntryPoint().y());
    return "\t" + record.join("  ") + "  /\n";
}
QString Welsegs::createSegmentEntry(Segment segment) {
 /*
  * 0. Segment number (first in range).
  * 1. Segment number (last in range).
  * 2. Branch number.
  * 3. Outlet segment number.
  * 4. Length of segment.
  * 5. Depth change through segment.
  * 6. Diameter.
  * 7. Roughness.
  * 8. * Cross sectional area for fluid flow.
  */
    auto record = GetBaseEntryLine(9);
    record[0] = QString::number(segment.Index());
    record[1] = QString::number(segment.Index());
    record[2] = QString::number(segment.Branch());
    record[3] = QString::number(segment.Outlet());
    record[4] = QString::number(segment.Length());
    record[5] = QString::number(segment.TVDChange());
    record[6] = QString::number(segment.Diameter());
    record[7] = QString::number(segment.Roughness());
    return "\t" + record.join("  ") + "  /";
}
Welsegs::Welsegs(QList<Model::Wells::Well *> *wells, int ts) {
    for (Well *well : *wells) {
        if (well->IsSegmented() && well->controls()->first()->time_step() == ts) {
            WelsegsKeyword kw;
            kw.heel_entry = createHeelEntry(well);
            auto segs = well->GetSegments();
            for (int i = 1; i < segs.size(); ++i) {
                kw.seg_entries.push_back(createSegmentEntry(segs[i]));
            }
            keywords_.push_back(kw);
        }
    }

}
QString Welsegs::WelsegsKeyword::buildKeyword() const {
    QString kw = "WELSEGS\n";
    kw += this->heel_entry;
    kw += this->seg_entries.join("\n");
    kw += "\n/\n\n";
    return kw;
}
}
}
