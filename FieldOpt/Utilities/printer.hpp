/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2017-2020 Mathias Bellout
<chakibbb-pcg@gmail.com>

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

// This file contains helper functions
// for "prettyfied" printing to console.

#ifndef PRINTER_FUNCTIONS_H
#define PRINTER_FUNCTIONS_H

#include <QString>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <stdio.h>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <Utilities/colors.hpp>
#include <Utilities/verbosity.h>

// ---------------------------------------------------------
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::to_string;
using std::stringstream;
using std::ofstream;
using std::setfill;
using std::left;
using std::setprecision;
using std::fixed;
using std::setw;
using std::runtime_error;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

namespace Printer {

inline void E(string m, string md, string cl) {
  m = "[mod: " + md + "] [cls: " + cl + "] " + m;
  throw runtime_error(m);
};

inline void DBG_prntToFile(string fn, string sout, const string mod="a") {
  FILE * pFile;
  pFile = std::fopen(fn.c_str(), mod.c_str());
  std::fprintf(pFile, "%s", sout.c_str());
  fclose (pFile);
}

inline string DBG_prntDbl(double argout,
                          const string& fd="% 10.3e ",
                          string fn="") {
  stringstream ss;
  char buffer [100];
  sprintf(buffer, fd.c_str(), argout);
  ss << buffer;
  // if (fn != "") DBG_printToFile(fn, ss.str() + "\n");
  return ss.str();
}

inline string DBG_prntVecXd(VectorXd vec, const string& mv="",
                            const string& fv="% 10.3e ",
                            string fn="", int sz=0) {
  stringstream ss;
  if ( vec.size() > 0 ) {
    if (mv != "" && mv != ",") ss << mv;
    ss << "[";
    for (int ii = 0; ii < vec.size() - 1; ii++) {
      ss << DBG_prntDbl(vec(ii), fv);
      if (mv == ",") ss << mv;
    }
    ss << DBG_prntDbl(vec(vec.size() - 1), fv) << "]";
    if (sz > 0) {
      ss << " [sz:" << DBG_prntDbl(vec.size(), "%2.0f") << "]";
    }
  } else {
    ss << "[ Empty ]";
  }

  // if (fn != "") DBG_prntToFile(fn, ss.str() + "\n");
  return ss.str();
}

inline string DBG_prntVecXdSz(VectorXd vec, const string& mv="",
                              const string& fv="% 10.3e ",
                              string fn="") {
  return DBG_prntVecXd(vec, mv, fv, fn, 1);
}

inline string DBG_prntMatXd(MatrixXd mat, string mm="",
                            string fm="% 10.3e ") {
  stringstream ss;
  for (int ii = 0; ii < mat.cols(); ii++) {
    if (mm != "") ss << "\n" << mm << "#" << ii;
    ss << DBG_prntVecXd(mat.col(ii), mm, fm);
  }
  return ss.str();
}

inline string RepStr(string str, int n) {
  std::ostringstream os;
  for(int i = 0; i < n; i++) { os << str; }
  return os.str();
}

struct BoxSym {

  string hll = "\u2500"; // ─ : horizontal line light
  string hlb = "\u2501"; // ━ : horizontal line bold
  string vll = "\u2502"; // ─ : vertical line light
  string vlb = "\u2503"; // ━ : vertical line bold

  string hdl = "\u2504"; // ─ : horizontal dashed-line light
  string hdb = "\u2505"; // ━ : horizontal dashed-line bold
  string vdl = "\u2506"; // ─ : vertical dashed-line light
  string vdb = "\u2507"; // ━ : vertical dashed-line bold

  string ull = "\u250C"; // ─ : upper left corner light
  string ulb = "\u250F"; // ━ : upper left corner bold
  string url = "\u2510"; // ─ : upper right corner light
  string urb = "\u2513"; // ━ : upper right corner bold

  string lll = "\u2514"; // ─ : lower left corner light
  string llb = "\u2517"; // ━ : lower left corner bold
  string lrl = "\u2518"; // ─ : lower right corner light
  string lrb = "\u251B"; // ━ : lower right corner bold

  string ulnl = ""; // upper line light
  string ulnb = ""; // upper line bold
  string llnl = ""; // lower line light
  string llnb = ""; // lower line bold

// Box Drawing symbols in unicode

//         0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F
// U+250x  ─   ━   │   ┃   ┄   ┅   ┆   ┇   ┈   ┉   ┊   ┋   ┌   ┍   ┎   ┏
// U+251x  ┐   ┑   ┒   ┓   └   ┕   ┖   ┗   ┘   ┙   ┚   ┛   ├   ┝   ┞   ┟
// U+252x  ┠   ┡   ┢   ┣   ┤   ┥   ┦   ┧   ┨   ┩   ┪   ┫   ┬   ┭   ┮   ┯
// U+253x  ┰   ┱   ┲   ┳   ┴   ┵   ┶   ┷   ┸   ┹   ┺   ┻   ┼   ┽   ┾   ┿
// U+254x  ╀   ╁   ╂   ╃   ╄   ╅   ╆   ╇   ╈   ╉   ╊   ╋   ╌   ╍   ╎   ╏
// U+255x  ═   ║   ╒   ╓   ╔   ╕   ╖   ╗   ╘   ╙   ╚   ╛   ╜   ╝   ╞   ╟
// U+256x  ╠   ╡   ╢   ╣   ╤   ╥   ╦   ╧   ╨   ╩   ╪   ╫   ╬   ╭   ╮   ╯
// U+257x  ╰   ╱   ╲   ╳   ╴   ╵   ╶   ╷   ╸   ╹   ╺   ╻   ╼   ╽   ╾   ╿

};

inline BoxSym bSym(int lw=60, bool dbg=false) {
  std::stringstream ss;
  BoxSym bS;

  bS.ulnl = bS.ull + RepStr(bS.hll, lw) + bS.url + "\n"; // upper line light
  bS.ulnb = bS.ull + RepStr(bS.hlb, lw) + bS.urb + "\n"; // upper line bold
  bS.llnl = bS.lll + RepStr(bS.hll, lw) + bS.lrl + "\n"; // lower line light
  bS.llnb = bS.lll + RepStr(bS.hlb, lw) + bS.lrb + "\n"; // lower line bold

  if(dbg) {
    ss << "horizontal line light:        " << bS.hll << " -- ";
    ss << "horizontal line bold:         " << bS.hlb << " -- ";
    ss << "vertical line light:          " << bS.vll << " -- ";
    ss << "vertical line bold:           " << bS.vlb << " -- \n";

    ss << "horizontal dashed-line light: " << bS.hdl << " -- ";
    ss << "horizontal dashed-line bold:  " << bS.hdb << " -- ";
    ss << "vertical dashed-line light:   " << bS.vdl << " -- ";
    ss << "vertical dashed-line bold:    " << bS.vdb << " -- \n";

    ss << "upper left corner light:      " << bS.ull << " -- ";
    ss << "upper left corner bold:       " << bS.ulb << " -- ";
    ss << "upper right corner light:     " << bS.url << " -- ";
    ss << "upper right corner bold:      " << bS.urb << " -- \n";

    ss << "lower left corner light:      " << bS.lll << " -- ";
    ss << "lower left corner bold:       " << bS.llb << " -- ";
    ss << "lower right corner light:     " << bS.lrl << " -- ";
    ss << "lower right corner bold:      " << bS.lrb << " -- \n";

    ss << "upper line light:             " << bS.ulnl << "\n";
    ss << "upper line bold:              " << bS.ulnb << "\n";
    ss << "lower line light:             " << bS.llnl << "\n";
    ss << "lower line bold:              " << bS.llnb << "\n\n";

    cout << ss.str();
  }

  return bS;
}

/*!
 * @brief Convert a number to its string representation.
 * @param num Number to convert.
 * @return String representation of num.
 */
template<typename T>
inline std::string num2str(const T num, int prc=2, int sci=0, int wdt=0) {
  // return boost::lexical_cast<std::string>(num);
  std::ostringstream ss;
  if (sci > 0) {
    if (wdt > 0) {
      ss << std::setw (wdt) << std::scientific << std::setprecision(prc) << num;
    } else {
      ss << std::scientific << std::setprecision(prc) << num;
    }
  } else {
    if (wdt > 0) {
      ss << std::setw (wdt) << std::fixed << std::setprecision(prc) << num;
    } else {
      ss << std::fixed << std::setprecision(prc) << num;
    }
  }
  return ss.str();
}

template<typename T>
inline QString num2strQ(const T num, int prc=2, int sci=0, int wdt=0) {
  return QString::fromStdString(num2str(num, prc, sci, wdt));
}

//template<typename T>
//inline QString num2Qstr(const T num, int prc=2) {
//  return QString::fromStdString(num2str(num, prc));
//}

/*!
 * @brief Truncate a string, adding an ellipsis to the end.
 * @param text Text to truncate.
 * @param width Width to to truncate to (including the added ellipsis).
 */
inline void truncate_text(string &text, const int &width=66) {
  if (text.size() > width) {
    text = text.substr(0, width - 3) + "...";
  }
}

/*!
 * @brief Pad a string, adding spaces at the
 * _end_ to make it the required length.
 * @param text Text to pad.
 * @param width Width to pad to.
 */
inline void pad_text(string &text, const int &width=66, const char &ps=' ') {
  if (text.size() < width) {
    text.insert(text.size(), width - text.size(), ps);
  }
}

/*!
 * @brief Split a string into multiple lines of the required
 * width, padding each with spaces at the end. Note that the
 * | character in the input string will insert a linebreak.
 * @param text Text to split.
 * @param width Width to split at (and pad to).
 * @return
 */
inline vector<string> split_line(const string text,
                                 const int &width=67) {
  vector<string> strv;
  if (text.size() > width) {
    std::size_t start_idx = 0;
    std::size_t end_idx = 1;
    int line_nr = 0;
    string remainder = text;
    while (remainder.size() > width
    || remainder.find_first_of("|") != std::string::npos) {
      if (remainder.find_first_of("|") < width) {
        end_idx = remainder.find_first_of("|");
      } else {
        end_idx = remainder.find_last_of(".,;/ ", width);
      }
      string line = remainder.substr(start_idx, end_idx);
      pad_text(line, width);
      strv.push_back(line);
      remainder = remainder.substr(end_idx+1, remainder.size());
      if (end_idx >= width) {
        break;
      }
    }
    pad_text(remainder, width);
    strv.push_back(remainder);
  } else {
    string padded = text;
    pad_text(padded, width);
    strv.push_back(padded);
  }
  return strv;
}


/*! @brief Print a compact infobox.
 Example:
 ┌─────────────────────────────────────────────────────────────────────┐
 │ This is a compact info box                                         │
 └─────────────────────────────────────────────────────────────────────┘
 */
inline void info(const std::string &text, int lw=163) {
  std::stringstream ss;
  ss << FLGREEN;
  std::string content = text;
  truncate_text(content, lw - 10);
  pad_text(content, lw - 4);

  BoxSym bS = bSym(lw);
  // ss << "\n" << bS.ulnl;
  ss << "\n";
  ss << "│ ■ " << content << " │" << "\n";
  // ss << bS.llnl;
  // ss << "\n";
  ss << AEND;
  std::cout << ss.str();
}

// ├───────┴──────────────────────┴──────────────────────────────────────┤
inline void idbg(const std::string &text, int lw=163) {
  std::stringstream ss;
  // ss << "\n"; // << FLGREEN
  auto lines = split_line(text, lw - 4);
  for (auto line : lines) {
    pad_text(line, lw - 4);
    // ss << "│ ■ " << line << " │" << "\n";
    ss << "\n" << "| = " << line << " |";
  }
  // ss << AEND;
  std::cout << ss.str();
}

/* Extended info box.
Example:
 ┌───────┬──────────────────────┬──────────────────────────────────────┐
 │ INFO │ Module: Optimization │ Class: Optimizer                     │
 ├───────┴──────────────────────┴──────────────────────────────────────┤
 │ This box can contain more information and a significantly larger    │
 │ amount of details.                                                  │
 └─────────────────────────────────────────────────────────────────────┘
 */
inline void ext_info(const std::string &text,
                     const std::string &modulen="",
                     const std::string &classn="",
                     const int &lw=163,
                     const int &lb=2,
                     const int &cl=0) {
  std::string module_name = modulen;
  std::string class_name = classn;
  truncate_text(module_name, 30);
  truncate_text(class_name, 30);
  pad_text(module_name, 30);
  pad_text(class_name, 30);

  auto lines = split_line(text, lw - lb);
  std::stringstream ss;

  if (cl==0) {
    ss << FLGREEN;
  } else if (cl==1) {
    ss << FLRED;
  } else if (cl==2) {
    ss << FLBLUE;
  } else if (cl==3) {
    ss << FLMAGENTA;
  }

  // auto width = floor(lw * 25/60);
  auto width = floor((lw-31) * 1/2) + 1;
  BoxSym bS = bSym(lw);
  ss << "\n" << bS.ulnl;
  pad_text(module_name, width);
  if (lw == 163) {
    pad_text(class_name, width);
  } else if (lw == 144) {
    pad_text(class_name, width+1);
  } else {
    pad_text(class_name, width+1);
  }

  ss << "│ ■ INFO │ Module: " << module_name << " │ Class: " << class_name << " │" << "\n";
  for (auto line : lines) {
    pad_text(line, lw - 2);
    ss << "│ " << line << " │" << "\n";
  }
  ss << bS.llnl;
  ss << "\n";
  ss << AEND;
  std::cout << ss.str();
}

/* Extended warning box.
Example:
 ┏━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
 ┃ WARNING  ┃ Module: Optimization ┃ Class: GeneticAlgorit...         ┃
 ┣━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┫
 ┃ This is a warning box.                                              ┃
 ┃ This is a warning box.                                              ┃
 ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
 */
inline void ext_warn(const std::string &text,
                     const std::string &modulen="",
                     const std::string &classn="",
                     const int &lw=163,
                     const int &cl=0) {
  std::string module_name = modulen;
  std::string class_name = classn;
  truncate_text(module_name, 30);
  truncate_text(class_name, 30);
  pad_text(module_name, 30);
  pad_text(class_name, 30);

  auto lines = split_line(text, lw - 2);
  std::stringstream ss;

  if (cl==0) {
    ss << FLYELLOW;
  } else if (cl==1) {
    ss << FLRED;
  } else if (cl==2) {
    ss << FLBLUE;
  } else if (cl==3) {
    ss << FLMAGENTA;
  }

  // auto width = floor(lw * 25/60);
  auto width = floor((lw-31) * 1/2) + 1;
  BoxSym bS = bSym(lw);
  ss << bS.ulnl;
  pad_text(module_name, width);
  if (lw == 163) {
    pad_text(class_name, width);
  } else if (lw == 144) {
    pad_text(class_name, width+1);
  } else {
    pad_text(class_name, width+1);
  }

  ss << "│ ■ WARN │ Module: " << module_name << " │ Class: " << class_name << " │" << "\n"; // #chars: 31
  for (auto line : lines) {
    pad_text(line, lw - 2);
    ss << "│ " << line << " │" << "\n";
  }
  ss << bS.llnl;
  ss << "\n";
  ss << AEND;
  std::cout << ss.str();
}

/* Error box.
Example:
 ╔═════════════════════════════════════════════════════════════════════╗
 ║ ERROR: This is an error message.                                   ║
 ╚═════════════════════════════════════════════════════════════════════╝
 */
inline void error(const std::string &text, const int &lw=163) {
  std::string content = text;

  auto lines = split_line(content, lw - 11);
  std::stringstream ss;
  ss << FLRED;

  BoxSym bS = bSym(lw);
  ss << bS.ulnl;
  for (auto line : lines) {
    pad_text(line, lw - 11);
    ss << "│ ■ ERROR: " << line << " │" << "\n";
  }
  ss << bS.llnl;
  ss << "\n";
  ss << AEND;
  std::cout << ss.str();
}

inline void ext_error(const std::string &text,
                      const std::string &modulen="",
                      const std::string &classn="",
                      const int &lw=163) {
  std::string module_name = modulen;
  std::string class_name = classn;
  truncate_text(module_name, 30);
  truncate_text(class_name, 30);
  pad_text(module_name, 30);
  pad_text(class_name, 30);

  auto lines = split_line(text, lw - 2);
  std::stringstream ss;
  ss << FRED;

  // auto width = floor(lw * 25/60);
  auto width = floor((lw-31) * 1/2) + 1;
  BoxSym bS = bSym(lw);
  ss << bS.ulnl;
  pad_text(module_name, width);
  pad_text(class_name, width+1);
  ss << "│ ■ WARN │ Module: " << module_name << " │ Class: " << class_name << " │" << "\n"; // #chars: 31
  for (auto line : lines) {
    pad_text(line, lw - 2);
    ss << "│ " << line << " │" << "\n";
  }
  ss << bS.llnl;
  ss << "\n";
  ss << AEND;
  std::cout << ss.str();
}

}


#endif // PRINTER_FUNCTIONS_H

