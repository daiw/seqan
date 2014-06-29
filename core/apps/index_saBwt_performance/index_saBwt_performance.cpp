// ==========================================================================
//                          index_saBwt_performance
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: David Iwanowitsch <iwanowit@inf.fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

#include <seqan/arg_parse.h>
#include "index_saBwt_performance.h"

// ==========================================================================
// Classes
// ==========================================================================



// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(SaBwtAppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("index_saBwt_performance");
    // Set short description, version, and date.
    setShortDescription(parser, "Measure Performance for index_sa_bpr and index_bwt_bcr");
    setVersion(parser, "0.1");
    setDate(parser, "July 2012");

    // Define usage line and long description.
    addUsageLine(parser, "[\\fIOPTIONS\\fP] \"\\fITEXT\\fP\"");
    addDescription(parser, "This is the application skelleton and you should modify this string.");

    // We require one argument.
    addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));
    addOption(parser, seqan::ArgParseOption("sa", "computeSA", "Compute Suffix Array"));
    addOption(parser, seqan::ArgParseOption("bwt", "computeBWT", "Compute Burrows-Wheeler-Transformation"));
    addOption(parser, seqan::ArgParseOption("comp", "compare", "Compare computation with refeerence computation"));
    addOption(parser, seqan::ArgParseOption("set", "stringSet", "InputFile contains a StringSet"));
    addOption(parser, seqan::ArgParseOption("f", "file", "Imputfile"));
    addOption(parser, seqan::ArgParseOption("ascii", "ascii", "input ist ascii not dna"));

    // Add Examples Section.
    addTextSection(parser, "Examples");
    addListItem(parser, "\\fBindex_saBwt_performance\\fP \\fB-v\\fP \\fItext\\fP",
                "Call with \\fITEXT\\fP set to \"text\" with verbose output.");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
    if (isSet(parser, "computeSA"))
       options.computeSA = true;
	if (isSet(parser, "computeBWT"))
		options.computeBwt = true;
	if (isSet(parser, "compare"))
		options.compare = true;
	if (isSet(parser, "set"))
		options.set = true;
	if (isSet(parser, "ascii"))
	    options.ascii = true;

	seqan::getArgumentValue(options.file, parser, 0);

    return seqan::ArgumentParser::PARSE_OK;
}

// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    SaBwtAppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "Params: "<<(options.computeSA?"Compute SA, ":"")<< (options.computeBwt?"Compute BWT, ":"")<<(options.compare?"Compare with reference, ":"")<<"File: "<<options.file<<std::endl;

    if(options.set){
    	doTheWorkSet(options);
    }else{
    	doTheWork(options);
    }

//    // Print the command line arguments back to the user.
//    if (options.verbosity > 0)
//    {
//        std::cout << "__OPTIONS____________________________________________________________________\n"
//                  << '\n'
//                  << "VERBOSITY\t" << options.verbosity << '\n'
//                  << "TEXT     \t" << options.text << "\n\n";
//    }
    std::cout<<std::endl;
    return 0;
}
