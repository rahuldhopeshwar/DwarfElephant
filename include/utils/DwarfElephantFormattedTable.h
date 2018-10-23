#ifndef DWARFELEPHANTFORMATTEDTABLE_H
#define DWARFELEPHANTFORMATTEDTABLE_H

#include "FormattedTable.h"

///Currently not used
class DwarfElephantFormattedTable : public FormattedTable
{
public:
  DwarfElephantFormattedTable();

  void printDakotaOutput(const std::string & file_name, int interval = 1, bool align = false);

private:
  std::ofstream _output_file;
  std::size_t _output_row_index;
  bool _stream_open;
  bool _append;
  std::string _csv_delimiter;
  unsigned int _csv_precision;
  std::string _output_file_name;

  void printRow(std::pair<Real, std::map<std::string, Real>> & row_data, bool align);
  void open(const std::string & file_name);
  void close();

  friend void
  dataStore<DwarfElephantFormattedTable>(std::ostream & stream, DwarfElephantFormattedTable & table, void * context);
  friend void dataLoad<DwarfElephantFormattedTable>(std::istream & stream, DwarfElephantFormattedTable & v, void * context);
};

template <>
void dataStore(std::ostream & stream, DwarfElephantFormattedTable & table, void * context);
template <>
void dataLoad(std::istream & stream, DwarfElephantFormattedTable & v, void * context);
#endif // DwarfElephantFORMATTEDTABLE_H
