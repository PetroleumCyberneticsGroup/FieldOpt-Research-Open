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

#include "compsegs.h"

namespace Simulation {
namespace ECLDriverParts {

using Printer::num2strQ;

Compsegs::Compsegs(Well *well) {
  head_ = "COMPSEGS\n";
  foot_ = "/\n\n";
  auto asegs = well->GetAnnulusSegments();
  CompsegsKeyword kw;
  kw.wname = well->name();
  for (int i = 0; i < asegs.size(); ++i) {
    kw.entries.push_back(generateEntry(asegs[i]));
  }
  keywords_.push_back(kw);
}

Compsegs::Compsegs(QList<Model::Wells::Well *> *wells, int ts) {
  for (Well *well : *wells) {
    if (well->IsSegmented() && well->controls()->first()->time_step() == ts) {
      CompsegsKeyword kw;
      kw.wname = well->name();

      auto segs = well->GetSegments();
      vector<Segment> tsegs;
      vector<Segment> isegs;
      vector<vector<Segment>> asegs(well->GetNumCompartments());

      for (const auto& seg : segs) {
        if (seg.Type() == Segment::TUBING_SEGMENT) {
          tsegs.push_back(seg);
        } else if (seg.Type() == Segment::ICD_SEGMENT) {
          isegs.push_back(seg);
        } else if (seg.Type() == Segment::ANNULUS_SEGMENT) {
          int k = seg.BlockBranch() - well->GetNumCompartments() - 2;
          asegs[k].push_back(seg);
        }
      }

      // dbg
      // cout << "tsegs.sz(): " << tsegs.size() << endl;
      // cout << "isegs.sz(): " << isegs.size() << endl;
      // cout << "asegs.sz(): " << asegs.size() << endl;
      // for (int kk=0; kk < asegs.size(); kk++) {
      //   cout << "asegs[kk=" << kk << "].sz(): " << asegs[kk].size() << endl;
      // }

      int idx_icd = 0, idx_ann = 0;

      for (int idx_tub=1; idx_tub < tsegs.size(); idx_tub++) {
        kw.entries.push_back(sepLine(idx_tub+1));
        kw.entries.push_back(generateEntry(tsegs[idx_tub]));

        kw.entries.push_back(generateEntry(isegs[idx_icd]));
        for (const auto & aseg : asegs[idx_icd]) {
          kw.entries.push_back(generateEntry(aseg));
        }
        idx_icd++;

      }

      // // ORIG CODE ----
      // // auto segs = well->GetSegments();
      // for (int i = 1; i < segs.size(); ++i) {
      //   kw.entries.push_back(generateEntry(segs[i]));
      // }
      //
      // // Some sort of sorting
      // QStringList qsl;
      // auto idx = getCompsegOrder(idx_orig_);
      // for (int i = 0; i < idx.size(); ++i) {
      //   // cout << "IDX=" << idx[i] << endl;
      //   // cout << kw.entries[idx[i]].toStdString() << endl;
      //   qsl << kw.entries[idx[i]];
      // }
      // // qSort(qsl.begin(), qsl.end());
      // // kw.entries = qsl;

      // ----
      // sorting but not by branch
      // qSort(kw.entries.begin(), kw.entries.end());
      keywords_.push_back(kw);

    }
  }
}

QString Simulation::ECLDriverParts::Compsegs::GetPartString() const {
  if (keywords_.size() == 0)
    return "";

  QString all_keywords = "";
  for (auto kw : keywords_) {
    all_keywords += kw.buildKeyword();
  }
  return all_keywords;
}

vector<size_t> Simulation::ECLDriverParts::Compsegs::getCompsegOrder(const vector<size_t> &v) {

  //  vector<size_t> idx(v.size());
  //  iota(idx.begin(), idx.end(), 0);
  // stable_sort(idx.begin(), idx.end(),
  //              [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);
  stable_sort(idx.begin(), idx.end(),
              [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

QString Compsegs::generateEntry(Segment seg) {
  /*
 * Record 1:
 *  1. Name of well.
 *
 * Remaining records:
 *  0. I-location
 *  1. J-location
 *  2. K-location
 *  3. Branch number.
 *  4. MD from root to _start_ of connection in this block.
 *  5. MD from root to _end_ of connection in this block.
 */
  auto entry = GetBaseEntryLine(6);
  // cout << "-> " << "I" << endl;
  entry[0] = num2strQ(seg.ParentBlock()->i(), 0, 0, 4) + "";
  // cout << "-> " << "J" << endl;
  entry[1] = num2strQ(seg.ParentBlock()->j(), 0, 0, 3) + "";
  // cout << "-> " << "K" << endl;
  entry[2] = num2strQ(seg.ParentBlock()->k(), 0, 0, 3) + "";
  // cout << "-> " << "Branch" << endl;
  entry[3] = num2strQ(seg.BlockBranch(), 0, 0, 4) + "   ";

  entry[4] = num2strQ(seg.GetEntryMD(), 5, 0, 9) + "";
  entry[5] = num2strQ(seg.GetExitMD(), 5, 0, 9) + "";

  // entry[4] = num2strQ(seg.OutletMD(), 5, 0, 9) + "";
  // entry[5] = num2strQ(seg.OutletMD() + seg.Length(), 5, 0, 9) + "";

  // entry[5] = num2strQ(seg.OutletMD() + seg.Length(), 5, 0, 9) + "  ";

  idx_orig_.push_back(seg.BlockBranch());
  return "" + entry.join("  ") + "  /";
}

QString Compsegs::CompsegsKeyword::buildKeyword() const {
  QString kw = "COMPSEGS\n";
  kw += "  " + this->wname + " /\n";
  kw += "-- I    J    K    Branch  Start     End        Dir Pen     End Range     Connection Depth\n";
  kw += "--                no      Length    Length     -\n";
  kw += this->entries.join("\n");
  kw += "\n/\n\n";
  return kw;
}
}
}
