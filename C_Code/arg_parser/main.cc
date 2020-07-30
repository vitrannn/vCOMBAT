/*  Arg_parser - POSIX/GNU command line argument parser. (C++ version)
    Copyright (C) 2006-2015 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
    Exit status: 0 for a normal exit, 1 for environmental problems
    (file not found, invalid flags, I/O errors, etc), 2 to indicate a
    corrupt or invalid input file, 3 for an internal consistency error
    (eg, bug) which caused arg_parser to panic.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "arg_parser.h"


namespace {

const char * const Program_name = "Arg_parser";
const char * const program_name = "arg_parser";
const char * const program_year = "2015";
const char * invocation_name = 0;


void show_help( const bool verbose )
  {
  std::printf( "%s - POSIX/GNU command line argument parser. (C++ version)\n", Program_name );
  std::printf( "See the source file 'main.cc' to learn how to use %s in\n", Program_name );
  std::printf( "your own programs.\n"
               "\nUsage: %s [options]\n", invocation_name );
  std::printf( "\nOptions:\n"
               "  -h, --help                   display this help and exit\n"
               "  -V, --version                output version information and exit\n"
               "  -a, --append                 example of option with no argument\n"
               "  -b, --block=<arg>            example of option with required argument\n"
               "  -c, --casual[=<arg>]         example of option with optional argument\n"
               "  -o <arg>                     example of short only option\n"
               "      --orphan                 example of long only option\n"
               "  -q, --quiet                  quiet operation\n"
               "  -u, --uncaught               example of intentional bug\n"
               "  -v, --verbose                verbose operation\n" );
  if( verbose )
    {
    std::printf( "  -H, --hidden                 example of hidden option (shown with -v -h)\n" );
    }
  std::printf( "\nReport bugs to arg-parser-bug@nongnu.org\n"
               "Arg_parser home page: http://www.nongnu.org/arg-parser/arg_parser.html\n" );
  }


void show_version()
  {
  std::printf( "%s %s\n", program_name, PROGVERSION );
  std::printf( "Copyright (C) %s Antonio Diaz Diaz.\n", program_year );
  std::printf( "License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n"
               "This is free software: you are free to change and redistribute it.\n"
               "There is NO WARRANTY, to the extent permitted by law.\n" );
  }


void show_error( const char * const msg, const int errcode = 0,
                 const bool help = false )
  {
  if( msg && msg[0] )
    {
    std::fprintf( stderr, "%s: %s", program_name, msg );
    if( errcode > 0 )
      std::fprintf( stderr, ": %s.", std::strerror( errcode ) );
    std::fprintf( stderr, "\n" );
    }
  if( help )
    std::fprintf( stderr, "Try '%s --help' for more information.\n",
                  invocation_name );
  }


void internal_error( const char * const msg )
  {
  std::fprintf( stderr, "%s: internal error: %s\n", program_name, msg );
  std::exit( 3 );
  }


const char * optname( const int code, const Arg_parser::Option options[] )
  {
  static char buf[2] = "?";

  if( code != 0 )
    for( int i = 0; options[i].code; ++i )
      if( code == options[i].code )
        { if( options[i].name ) return options[i].name; else break; }
  if( code > 0 && code < 256 ) buf[0] = code; else buf[0] = '?';
  return buf;
  }

} // end namespace


int main( const int argc, const char * const argv[] )
  {
  bool verbose = false;
  invocation_name = argv[0];

  const Arg_parser::Option options[] =
    {
    { 'a', "append",   Arg_parser::no    },
    { 'b', "block",    Arg_parser::yes   },
    { 'c', "casual",   Arg_parser::maybe },
    { 'h', "help",     Arg_parser::no    },
    { 'H', "hidden",   Arg_parser::no    },
    { 'o', 0,          Arg_parser::yes   },
    { 'q', "quiet",    Arg_parser::no    },
    { 'u', "uncaught", Arg_parser::no    },
    { 'v', "verbose",  Arg_parser::no    },
    { 'V', "version",  Arg_parser::no    },
    { 256, "orphan",   Arg_parser::no    },
    {   0, 0,          Arg_parser::no    } };

  const Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  for( int argind = 0; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;				// no more options
    switch( code )
      {
      case 'a': break;				// example, do nothing
      case 'b': break;				// example, do nothing
      case 'c': break;				// example, do nothing
      case 'h': show_help( verbose ); return 0;
      case 'H': break;				// example, do nothing
      case 'o': break;				// example, do nothing
      case 'q': verbose = false; break;
      // case 'u': break;			// intentionally not caught
      case 'v': verbose = true; break;
      case 'V': show_version(); return 0;
      case 256: break;				// example, do nothing
      default : internal_error( "uncaught option." );
      }
    } // end process options

  for( int argind = 0; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    const char * const arg = parser.argument( argind ).c_str();
    if( code )	// option
      {
      const char * const name = optname( code, options );
      if( !name[1] )
        std::printf( "option '-%c'", name[0] );
      else
        std::printf( "option '--%s'", name );
      if( arg[0] )
        std::printf( " with argument '%s'", arg );
      }
    else	// non-option
      std::printf( "non-option argument '%s'", arg );
    std::printf( "\n" );
    }

  if( !parser.arguments() ) std::printf( "Hello, world!\n" );

  return 0;
  }
