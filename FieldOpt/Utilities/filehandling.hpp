/***********************************************************
Copyright (C) 2015-2017
Einar J.M. Baumann <einar.baumann@gmail.com>

Modified 2020 Mathias Bellout
<chakibbb.pcg@gmail.com>

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

#ifndef FILEHANDLING_H
#define FILEHANDLING_H

#include <QString>
#include <QStringList>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QDir>
#include <stdexcept>
#include <string>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <utility>
#include <vector>

#include "Utilities/verbosity.h"
#include <Utilities/printer.hpp>

namespace Utilities {
namespace FileHandling {

using Printer::ext_info;
using Printer::info;
using Printer::ext_warn;
using boost::filesystem::copy_file;

/*!
 * \brief FileExists Checks whether or not file exists at specified path.
 * \param file_path Path to a file that may or may not exist.
 * \param verbose Whether the path being checked should be printed.
 * \return True if a file exists at the specified path, otherwise false.
 */
inline bool FileExists(const QString &file_path, Settings::VerbParams vp,
                       string mod="", string cls="") {
  QFileInfo file(file_path);
  QFileInfo file_relative(file.absoluteFilePath());
  auto fp = file_path.toStdString();
  if (file.exists() && file.isFile()) {
    if (vp.vSIM >= 5) { ext_info("File exists at path: " + fp, mod, cls, vp.lnw); }
    return true;

  } else if (file_relative.exists() && file_relative.isFile()) {
    if (vp.vSIM >= 5) { ext_info("File exists at rel.path: " + fp, mod, cls, vp.lnw); }
    return true;

  } else {
    if (vp.vSIM >= 5) { ext_info("File does not exists: " + fp, mod, cls, vp.lnw); }
    return false;
  }
}

/*!
 * Overload of the FileExists(QString, bool) function. Creates a
 * QString from the std string and calls the other function.
 */
inline bool FileExists(const std::string& file_path,
                       const Settings::VerbParams vp,
                       string mod="", string cls="") {
  return FileExists(QString::fromStdString(file_path), vp,
                    std::move(mod), std::move(cls));
}

/*!
 * \brief DirExists Checks whether or not a folder exists at the specified path.
 * \param folder_path Path to a folder that may or may not exist.
 * \param verbose Whether the path being checked should be printed.
 * \return True if a folder exists at the specified path, otherwise false.
 */
inline bool DirExists(const QString &directory_path,
                      const Settings::VerbParams vp,
                      string mod="", string cls="") {
  QFileInfo folder(directory_path);
  auto dp = directory_path.toStdString();
  if (folder.exists() && folder.isDir()) {
    if (vp.vSIM >= 5) {
      auto tm = "Dir exists at path: " + dp;
      ext_info(tm, mod, cls, vp.lnw); }
    return true;

  } else {
    if (vp.vSIM >= 5) {
      auto tm = "Dir does not exists at path: " + dp;
      ext_info(tm, mod, cls, vp.lnw); }
    return false;
  }
}

/*!
 * Overload of the DirExists(QString, bool) function. Creates a
 * QString from the std string and calls the other function.
 */
inline bool DirExists(const std::string &directory_path,
                      const Settings::VerbParams vp,
                      string mod="", string cls="") {
  return DirExists(QString::fromStdString(directory_path),
                   vp, std::move(mod), std::move(cls));
}


/*!
 * \brief DirIsEmpty Check whether or not a directory is empty.
 * \param folder_path Path a folder to check.
 * \return True if the directory is empty, otherwise false.
 */
inline bool DirIsEmpty(const QString &folder_path,
                       const Settings::VerbParams vp) {
  if (!DirExists(folder_path, vp)) return false;
  QDir directory = QDir(folder_path);
  return directory.entryInfoList(QDir::NoDotAndDotDot |
      QDir::AllEntries).count() == 0;
}

inline bool DirIsEmpty(const std::string &folder_path,
                       const Settings::VerbParams vp) {
  return DirIsEmpty(QString::fromStdString(folder_path), vp);
}

/*!
 * \brief ParentDirExists Checks whether a specified file's parent directory exists.
 * \param file_path Path a file (the file itself does not have to exist).
 * \return True if the parent directory exists, otherwise false.
 */
inline bool ParentDirExists(QString file_path) {
  QFileInfo file(file_path);
  QDir parent_directory = file.dir();
  return parent_directory.exists();
}

/*!
 * \brief ReadFileToStringList Reads the contents of a file and stores it as
 * a string list where each element is a line in the file.
 * \param file_path The file to create a list from.
 * \return List where each element is a line in the file.
 */
inline QStringList *ReadFileToStringList(const QString &file_path) {
  QStringList *string_list = new QStringList();
  QFile file(file_path);
  if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    throw std::runtime_error("File not found: " + file_path.toStdString());

  QTextStream text_stream(&file);
  while (true) {
    QString line = text_stream.readLine();
    if (line.isNull())
      break;
    else
      string_list->append(line);
  }
  file.close();
  return string_list;
}

inline std::vector<std::string> ReadFileToStdStringList(const std::string &filepath) {
  auto qt_string_list = ReadFileToStringList(QString::fromStdString(filepath));
  auto stringlist = std::vector<std::string>(qt_string_list->size());

  for (int i = 0; i < qt_string_list->size(); ++i) {
    stringlist[i] = qt_string_list->at(i).toStdString();
  }
  qt_string_list->clear();
  delete qt_string_list;

  return stringlist;
}

/*!
 * \brief WriteStringToFile Write a string to a file. Removes existing file contents.
 *
 * If the string does not end with a newline, it will be added.
 * \param string The string to be written.
 * \param file_path Path to the file to write the string into.
 */
inline void WriteStringToFile(QString string, QString file_path) {
  std::string em;
  if (!ParentDirExists(file_path)) {
    em = "File's parent directory not found: " + file_path.toStdString();
    throw std::runtime_error(em);
  }

  if (!string.endsWith("\n"))
    string.append("\n");

  QFile file(file_path);
  file.open(QIODevice::WriteOnly | QIODevice::Truncate);
  QTextStream out(&file);
  out << string.toUtf8() << endl;
  file.close();
}

inline void WriteStringToFile(string str, string file_path) {
  WriteStringToFile(QString::fromStdString(str),
                    QString::fromStdString(file_path));
}

inline void WriteStringToFile(QString str, string file_path) {
  WriteStringToFile(str, QString::fromStdString(file_path));
}

/*!
 * \brief WriteLineToFile Append a string to a file.
 *
 * If the string does not end with a newline, it will be added.
 * \param string The string/line to be written.
 * \param file_path The file to write the string/line to.
 */
inline void WriteLineToFile(QString string, const QString& file_path) {
  if (!ParentDirExists(file_path))
    throw std::runtime_error("File's parent directory not found: " + file_path.toStdString());

  if (!string.endsWith("\n"))
    string.append("\n");

  QFile file(file_path);
  file.open(QIODevice::Append);
  QTextStream out(&file);
  out << string.toUtf8();
  file.close();
}

/*!
 * \brief DeleteFile Deletes the file at the given path.
 * \param path Path to file to be deleted.
 */
inline void DeleteFile(const QString& path,
                       const Settings::VerbParams vp,
                       string mod="", string cls="") {
  if (FileExists(path, vp, mod, cls)) {
    QFile file(path);
    file.remove();
  } else {
    throw std::runtime_error("File not found: " + path.toStdString());
  }
}

/*!
 * \brief CreateDir Create a new drectory with the specified path.
 * \param path Path to new directory.
 */
inline void CreateDir(QString path,
                      const Settings::VerbParams vp,
                      string mod="", string cls="") {
  if (DirExists(path, vp, mod, cls))
    return; // Do nothing if the directory already exists.
  QDir().mkdir(path);
}

inline void CreateDir(std::string path,
                      const Settings::VerbParams vp,
                      string mod="", string cls="") {
  CreateDir(QString::fromStdString(path), vp, mod, cls);
}

/*!
 * Get the name of a file from a path (i.e. delete everyting up to
 * and including the final /).
 * @param file_path Path to a file
 * @return Name of a file, including extension.
 */
inline std::string FileName(const std::string file_path) {
  std::vector<std::string> parts;
  boost::split(parts, file_path,
               boost::is_any_of("/"), boost::token_compress_on);
  return parts.back();
}

inline std::string FileNameRoot(const std::string file_path) {
  std::vector<std::string> parts;
  boost::split(parts, file_path,
               boost::is_any_of("."), boost::token_compress_on);
  return parts.front();
}

inline QString FileNameQstr(const std::string file_path) {
  return QString::fromStdString(FileName(file_path));
}

/*!
 * Get the name of a file's parent directory.
 * @param file_path Path to a file.
 * @return Name of a directory.
 */
inline std::string ParentDirName(const std::string file_path) {
  std::vector<std::string> parts;
  boost::split(parts, file_path, boost::is_any_of("/"),
               boost::token_compress_on);
  return parts[parts.size() - 2];
}

inline QString ParentDirectoryNameQstr(const std::string file_path) {
  return QString::fromStdString(ParentDirName(file_path));
}

/*!
 * \brief CopyFile Copy a file.
 * \param origin The path to the original file.
 * \param destination Path to the copy of the file.
 * \param overwrite Overwrite existing file.
 */
inline void CopyFile(QString origin, QString destination,
                     bool overwrite, Settings::VerbParams vp,
                     string mod="", string cls="") {
  if (!FileExists(origin, vp, mod, cls)) {
    throw std::runtime_error("Error copying. Original file not found: " + origin.toStdString());
  }

  if (overwrite) {
    copy_file(origin.toStdString(),
              destination.toStdString(),
              boost::filesystem::copy_option::overwrite_if_exists);
  } else {
    copy_file(origin.toStdString(), destination.toStdString());
  }
}

inline void CopyFile(const std::string& origin,
                     const std::string& destination,
                     bool overwrite=false,
                     Settings::VerbParams vp={},
                     string mod="", string cls="") {
  CopyFile(QString::fromStdString(origin),
           QString::fromStdString(destination),
           overwrite, vp, mod, cls);
}

/*!
 * \brief CopyDir
 * Copy a directory and it's contents to a new destination.
 *
 * Note this is not a recursive function: It will copy files
 * in the root of the directory, and _create_ any subdirectories
 * found, but it will not copy the contents of subdirectories.
 * \param origin Path to the original directory to be copied.
 * \param dest Path to the _parent directory_ for the copy.
 */
inline void CopyDir(QString origin, QString dest,
                    bool verbose=false, bool ecl_slim=false,
                    const Settings::VerbParams vp={},
                    string mod="", string cls="") {

  if (!DirExists(origin, vp)) {
    throw std::runtime_error("Parent dir for copying not found:\n"
                                 + origin.toStdString());
  }

  if (!DirExists(dest, vp)) {
    throw std::runtime_error("Destination (parent) dir for copying not found:\n"
                                 + dest.toStdString());
  }

  QDir original(origin);
  QFileInfoList entries;
  entries = original.entryInfoList(QDir::AllEntries | QDir::NoDotAndDotDot, QDir::DirsLast);

  QStringList eclOutFiles;
  eclOutFiles << "CONF1" << "CONF2" << "DBG" << "ECLEND" << "EGRID"
              << "GRID" << "INIT" << "INSPEC" << "MSG" << "OPT"
              << "OUT" << "PRT" << "RSM" << "RSSPEC" << "SMRY"
              << "SMRY_DBG" << "SMSPEC" << "SOLN_RUNT" << "SOLN_SCSA"
              << "SOLN_WBHP" << "UNSMRY" << "UNRST";

  for (auto entry : entries) {
    if (entry.isFile() && !entry.isDir()) {
      if (ecl_slim) {
        if (! eclOutFiles.contains(entry.suffix())) {
          CopyFile(entry.absoluteFilePath(),dest + "/" + entry.fileName(), true, vp, mod, cls);
          if (verbose) std::cout << "Copying FILE: " << entry.fileName().toStdString() << std::endl;
        }
      } else {
        CopyFile(entry.absoluteFilePath(), dest + "/" + entry.fileName(), true, vp, mod, cls);
        if (verbose) std::cout << "Copying FILE: " << entry.fileName().toStdString() << std::endl;
      }

    } else if (entry.isDir()) {
      CreateDir(dest + "/" + entry.fileName(), vp);
      if(verbose) std::cout << "Copying FOLDER: " << entry.fileName().toStdString() << std::endl;
      CopyDir(entry.absoluteFilePath(), dest + "/" + entry.fileName(), verbose, true, vp, mod, cls);
    }
  }
}

inline void CopyDir(const std::string& origin,
                    const std::string& destination,
                    bool verbose= false,
                    bool ecl_slim= false,
                    const Settings::VerbParams vp={},
                    string mod="", string cls="") {
  CopyDir(QString::fromStdString(origin),
          QString::fromStdString(destination),
          verbose, ecl_slim, vp, mod, cls);
}

/*!
 * \brief GetCurrentDirectoryPath
 * Gets the absolute path to the current directory.
 * \todo Improve this.
 */
inline QString GetCurrentDirectoryPath() {
  QDir path = QDir::currentPath(); // Get current directory
  return path.absolutePath();
}

/*!
 * \brief GetAbsoluteFilePath Gets absolute path of file.
 * \param file (relative) path to file
 */
inline QString GetAbsoluteFilePath(const QString &file) {
  QFileInfo fileInfo(file);
  return fileInfo.absoluteFilePath();
}

inline std::string GetAbsoluteFilePath(const std::string &file) {
  return GetAbsoluteFilePath(QString::fromStdString(file)).toStdString();
}

/*!
 * Get the path to a file's parent directory (i.e. remove everyting
 * after the final slash)
 */
inline QString GetParentDirPath(const QString &file_path) {
  QStringList parts = file_path.split("/");
  parts.removeLast();
  return parts.join("/");
}

inline std::string GetParentDirPath(const std::string &file_path) {
  return GetParentDirPath(QString::fromStdString(file_path)).toStdString();
}

}
}

#endif // FILEHANDLING_H
