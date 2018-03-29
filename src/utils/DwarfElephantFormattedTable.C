#include "DwarfElephantFormattedTable.h"
#include "MooseError.h"
#include "InfixIterator.h"
#include "MooseUtils.h"

#include "libmesh/exodusII_io.h"

#include <iomanip>
#include <iterator>

// Used for terminal width
#include <sys/ioctl.h>
#include <cstdlib>

const unsigned short DEFAULT_CSV_PRECISION = 14;
const std::string DEFAULT_CSV_DELIMITER = "   ";

template <>
void
dataStore(std::ostream & stream, DwarfElephantFormattedTable & table, void * context)
{
  storeHelper(stream, table._data, context);
  storeHelper(stream, table._align_widths, context);
  storeHelper(stream, table._column_names, context);
  storeHelper(stream, table._output_row_index, context);
}

template <>
void
dataLoad(std::istream & stream, DwarfElephantFormattedTable & table, void * context)
{
  loadHelper(stream, table._data, context);
  loadHelper(stream, table._align_widths, context);
  loadHelper(stream, table._column_names, context);
  loadHelper(stream, table._output_row_index, context);

  // Don't assume that the stream is open if we've restored.
  table._stream_open = false;
}


DwarfElephantFormattedTable::DwarfElephantFormattedTable()
  : FormattedTable(),
   _output_row_index(0),
   _stream_open(false),
   _append(false),
  _csv_delimiter(DEFAULT_CSV_DELIMITER),
  _csv_precision(DEFAULT_CSV_PRECISION)
{
}

void
DwarfElephantFormattedTable::printDakotaOutput(const std::string & file_name, int interval, bool align)
{
  open(file_name);

  if (_output_row_index == 0)
  {
    /**
     * When the alignment option is set to true, the widths of the columns needs to be computed
     * based on longest of the column name of the data supplied. This is done here by creating a
     * map
     * of the widths for each of the columns, including time
     */
    if (align)
    {
      // Set the initial width to the names of the columns
      _align_widths["time"] = 4;

      for (const auto & col_name : _column_names)
        _align_widths[col_name] = col_name.size();

      // Loop through the various times
      for (const auto & it : _data)
      {
        // Update the time _align_width
        {
          std::ostringstream oss;
          oss << std::setprecision(_csv_precision) << it.first;
          unsigned int w = oss.str().size();
          _align_widths["time"] = std::max(_align_widths["time"], w);
        }

        // Loop through the data for the current time and update the _align_widths
        for (const auto & jt : it.second)
        {
          std::ostringstream oss;
          oss << std::setprecision(_csv_precision) << jt.second;
          unsigned int w = oss.str().size();
          _align_widths[jt.first] = std::max(_align_widths[jt.first], w);
        }
      }
    }

    // Output Header
    {
      for (const auto & col_name : _column_names)
      {
        _output_file << _csv_delimiter;

        if (align)
          _output_file << std::right << std::setw(_align_widths[col_name]) << col_name;
        else
          _output_file << col_name;
      }
      _output_file << "\n";
    }
  }

  for (; _output_row_index < _data.size(); ++_output_row_index)
  {
    if (_output_row_index % interval == 0)
      printRow(_data[_output_row_index], align);
  }

  _output_file.flush();
}


void
DwarfElephantFormattedTable::printRow(std::pair<Real, std::map<std::string, Real>> & row_data, bool align)
{
  for (const auto & col_name : _column_names)
  {
    std::map<std::string, Real> & tmp = row_data.second;

    _output_file << _csv_delimiter;

    if (align)
      _output_file << std::setprecision(_csv_precision) << std::right
                   << std::setw(_align_widths[col_name]) << tmp[col_name];
    else
      _output_file << std::setprecision(_csv_precision) << tmp[col_name];
  }
  _output_file << "\n";
}

void
DwarfElephantFormattedTable::open(const std::string & file_name)
{
  if (_stream_open && _output_file_name == file_name)
    return;
  close();
  _output_file_name = file_name;

  std::ios_base::openmode open_flags = std::ios::out;
  if (_append)
    open_flags |= std::ios::app;
  else
  {
    open_flags |= std::ios::trunc;
    _output_row_index = 0;
  }

  _output_file.open(file_name.c_str(), open_flags);
  _stream_open = true;
}

void
DwarfElephantFormattedTable::close()
{
  if (!_stream_open)
    return;
  _output_file.flush();
  _output_file.close();
  _stream_open = false;
  _output_file_name = "";
}
