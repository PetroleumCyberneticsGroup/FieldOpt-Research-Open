/***********************************************************
Copyright (C) 2015-2018
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020-2021 Mathias Bellout
<mathias.bellout@gmail.com>

This file is part of the FieldOpt project.

FieldOpt is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version
3 of the License, or (at your option) any later version.

FieldOpt is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
the GNU General Public License for more details.

You should have received a copy of the
GNU General Public License along with FieldOpt.
If not, see <http://www.gnu.org/licenses/>.
***********************************************************/

#include "welsegs.h"

namespace Simulation {
namespace ECLDriverParts {

using Printer::num2strQ;

Welsegs::Welsegs(Well *well) {
  head_ = "WELSEGS\n";
  foot_ = "/\n\n";

  WelsegsKeyword kw;
  kw.heel_entry = createHeelEntry(well);
  auto segs = well->GetSegments();

  for (int i = 1; i < segs.size(); ++i) {
    if (segs[i].Type() == Segment::TUBING_SEGMENT) {
      kw.seg_entries.push_back(createSegmentEntry(segs[i]));
    }
  }

  kw.seg_entries.push_back("  " + sepLine());

  for (int i = 1; i < segs.size(); ++i) {
    if (segs[i].Type() == Segment::ICD_SEGMENT) {
      kw.seg_entries.push_back(createSegEntryICD(segs[i]));
    }
  }

  // kw.seg_entries.push_back("  " + sepLine());

  for (int i = 1; i < segs.size(); ++i) {
    if (segs[i].Type() == Segment::ANNULUS_SEGMENT) {
      kw.seg_entries.push_back(createSegmentEntry(segs[i]));
    }
  }
  keywords_.push_back(kw);
}

Welsegs::Welsegs(QList<Model::Wells::Well *> *wells, int ts) {
  for (Well *well : *wells) {
    if (well->IsSegmented() && well->controls()->first()->time_step() == ts) {

      WelsegsKeyword kw;
      kw.heel_entry = createHeelEntry(well);
      auto segs = well->GetSegments();

      for (int i = 1; i < segs.size(); ++i) {
        if (segs[i].Type() == Segment::TUBING_SEGMENT) {
          kw.seg_entries.push_back(createSegmentEntry(segs[i]));
        }
      }

      kw.seg_entries.push_back("  " + sepLine());

      for (int i = 1; i < segs.size(); ++i) {
        if (segs[i].Type() == Segment::ICD_SEGMENT) {
          kw.seg_entries.push_back(createSegEntryICD(segs[i]));
        }
      }

      // kw.seg_entries.push_back("  " + sepLine());

      for (int i = 1; i < segs.size(); ++i) {
        if (segs[i].Type() == Segment::ANNULUS_SEGMENT) {
          if (segs[i].AuxIndex() > 0) {
            kw.seg_entries.push_back("  " + sepLine(segs[i].AuxIndex()));
          }
          kw.seg_entries.push_back(createSegmentEntry(segs[i]));
        }
      }
      keywords_.push_back(kw);
    }
  }
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
  auto record = GetBaseEntryLine(11);

  record[0]  = "-- Name   Dep 1         Tlen 1    Vol 1    Len&Dep    PresDrop\n";
  record[1]  = "   " + well->name() + "  ";
  record[2]  = num2strQ(rseg.TVDChange(), 5, 0, 9) + "    ";
  record[3]  = num2strQ(0.0, 5, 0, 7) + "    ";
  record[5]  = "      INC  ";
  record[6]  = "      HF-  /\n";
  record[7]  = "--                                                                          -\n";
  record[8]  = "-- 1st    Last    Branch Outlet Length      Depth     Diam      Rough       -\n";
  record[9]  = "-- Seg    Seg     Num    Seg                Change                          -\n";
  record[10] = "-- Main Stem Segments                                                       -\n";
  // record[6] = "DF";
  // record[7] = QString::number(well->trajectory()->GetWellBlocks()->at(0)->getEntryPoint().x());
  // record[8] = QString::number(well->trajectory()->GetWellBlocks()->at(0)->getEntryPoint().y());
  return record.join("");
}

QString Welsegs::createSegEntryICD(Segment segment) {
  auto record = GetBaseEntryLine(10);
  record[0] = num2strQ(segment.Index(), 0, 0, 4);
  record[1] = "  " + num2strQ(segment.Index(), 0, 0, 3);
  record[2] = "  " + num2strQ(segment.Branch(), 0, 0, 3);
  record[3] = "  " + num2strQ(segment.Outlet(), 0, 0, 2);

  record[4] = "    " + num2strQ(segment.Length(), 5, 0, 9);
  record[5] = " " + num2strQ(segment.TVDChange(), 5, 0, 7);
  record[6] = " " + num2strQ(segment.Diameter(), 5, 0, 7);
  record[7] = " " + num2strQ(segment.Roughness(), 3, 1, 7);
  return "  " + record.join("  ") + "  /";
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
  auto record = GetBaseEntryLine(8);
  record[0] = num2strQ(segment.Index(), 0, 0, 4);
  record[1] = "  " + num2strQ(segment.Index(), 0, 0, 3);
  record[2] = "  " + num2strQ(segment.Branch(), 0, 0, 3);
  record[3] = "  " + num2strQ(segment.Outlet(), 0, 0, 2);

  record[4] = "    " + num2strQ(segment.Length(), 5, 0, 9);
  record[5] = " " + num2strQ(segment.TVDChange(), 5, 0, 7);
  record[6] = " " + num2strQ(segment.Diameter(), 5, 0, 7);
  record[7] = " " + num2strQ(segment.Roughness(), 5, 0, 7);
  return "  " + record.join("  ") + "  /";
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
