#ifndef DRIVERPART_H
#define DRIVERPART_H

#include <QString>
#include <QStringList>

namespace Simulation {

/*!
 * \brief The DriverPart class Is the parent class for any part of a driver file.
 */
class DriverPart
{
public:
    /*!
     * \brief GetPartString Get a string containing the section of the driver file.
     * \return
     */
    virtual QString GetPartString() const = 0;

protected:
    DriverPart() {}

    /*!
     * \brief getKeywordContent Get the content of the section in the driver file, i.e. everything from (and including)
     * the keyword and up to (but not including) the next keyword.
     * \param driver_file_content The driver file to search through.
     * \param keyword The keyword to get the contents of.
     * \param next_keyword The next keyword after the one you want.
     * \return The matching lines, joined with a newline between them.
     */
    QString getSectionContent(QStringList *driver_file_content, QString keyword, QString next_keyword);

    /*!
     * \brief getKeywordContent Get the content of the section in the driver file, i.e. everything from (and including)
     * the keyword and up to (but not including) the next keyword. This method is used if there are more than one possible
     * next keywords, i.e. when the next keyword can be an optional one.
     * \param driver_file_content The driver file to search through.
     * \param keyword The keyword to get the contents of.
     * \param possible_next_keywords A list of possible next keywords.
     * \return
     */
    QString getSectionContent(QStringList *driver_file_content, QString keyword, QStringList possible_next_keywords);

private:
    /*!
     * \brief getLines Get the lines from (and including) one index to (but not including) another index in the
     * string list.
     * \param driver_file_content List to get lines from.
     * \param start_index Index of first line to get.
     * \param end_index Index of the last line to get.
     * \return
     */
    QString getLines(QStringList *driver_file_content, int start_index, int end_index);

    /*!
     * \brief getLineIndex Get the index of the line on which the input keyword first occurs.
     * \param driver_file_content The lines to search through.
     * \param keyword The keyword to find the index of.
     * \return The index of the first line on which the keywrod occurs. Returns -1 if it is not found.
     */
    int getLineIndex(QStringList *driver_file_content, QString keyword);
};

}

#endif // DRIVERPART_H
