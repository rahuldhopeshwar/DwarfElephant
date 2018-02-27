#include "DwarfElephantApp.h"
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("DwarfElephant");

// Begin the main program.
int main(int argc, char *argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  DwarfElephantApp::registerApps();

  // In case our are using a MOOSE version older than Jan 24, 2018
  // use the following line
  // MooseApp * app = AppFactory::createApp("DwarfElephantApp", argc, argv);

  // Create an instance of the application and store it in a smart pointer for easy cleanup
  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("DwarfElephantApp", argc, argv);

  // Execute the application
  app->run();

  // In case our are using a MOOSE version older than Jan 24, 2018
  // use the following line
  // delete app;

  return 0;
}
