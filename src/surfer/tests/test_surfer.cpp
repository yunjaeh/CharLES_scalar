#define CATCH_CONFIG_RUNNER
#include "../../core/catch1/catch.hpp"

int main( int argc, char* argv[] )
{
  // global setup...
  
  int result = Catch::Session().run( argc, argv );
  
  // global clean-up...

  return ( result < 0xff ? result : 0xff );
}

