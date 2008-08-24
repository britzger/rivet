#include "Rivet/Rivet.hh"
#include "Rivet/HistoFormat.hh"
#include "Rivet/AnalysisHandler.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Tools/Configuration.hh"

#include "AGILe/Particle.hh"
#include "AGILe/Generator.hh"
#include "AGILe/Loader.hh"

#include <tclap/CmdLine.h>
using namespace TCLAP;

namespace {
  using namespace std;
  void handleConfigStream(istream& in, Rivet::ParamMap& pmap);
  pair<string,string> parseParamString(const string& paramstring);
}

#ifndef RIVETVERSION
#define RIVETVERSION "UNKNOWN"
#endif

///////////////////////////////////////////


namespace Rivet {
  namespace Commandline {


    Configuration parseArgs(size_t argc, char** argv) {
      // The Configuration to be returned
      Configuration config;

      try {
        CmdLine cmd("Runs an event generator using the RivetGun and AGILe interface libraries", ' ', RIVETVERSION);

        // Add generator args
        vector<string> gens = AGILe::Loader::getAvailableGens();
        ValuesConstraint<string> genNameConstraint(gens);
        vector<Arg*> genArgs;
        ValueArg<string> genNameArg("g", "generator", "Generator to use", true, "", &genNameConstraint);
        genArgs.push_back(&genNameArg);
        ValueArg<string> genFileArg("G", "generatorfile", "HepML file defining the generator", true, "", "genfile");
        genArgs.push_back(&genFileArg);
        ValueArg<string> eventFileArg("i", "ineventfile", "File containing HepMC events", true, "-", "filename");
        genArgs.push_back(&eventFileArg);
        //cmd.xorAdd(genNameArg, genFileArg);
        cmd.xorAdd(genArgs);
        ValueArg<size_t> rngSeedArg("s", "seed", "random number generator seed", false, 271828, "seed", cmd);

        // Add initial state args
        //vector<string> particles = Rivet::getParticleNames();
        vector<string> particles;
        particles += "PROTON","ANTIPROTON","ELECTRON","POSITRON";
        ValuesConstraint<string> beamNameConstraint(particles);
        ValueArg<string> beam1Arg("", "beam1", "Particle in beam 1 (PROTON by default)", false, "PROTON", &beamNameConstraint, cmd);
        ValueArg<string> beam2Arg("", "beam2", "Particle in beam 2 (PROTON by default)", false, "PROTON", &beamNameConstraint, cmd);
        ValueArg<double> mom1Arg("", "mom1", "Momentum of beam 1 in GeV (7000 by default)", false, 7000.0, "mom", cmd);
        ValueArg<double> mom2Arg("", "mom2", "Momentum of beam 2 in GeV (7000 by default)", false, 7000.0, "mom", cmd);

        // Add generator param args
        MultiArg<string> paramsArg("p", "param", "Param in 'name:value' format", false, "param=value", cmd);
        MultiArg<string> paramsFilesArg("P", "paramsfile", "File containing parameter name, value pairs", false, "paramfile", cmd);

        // Add analysis args
        const set<string> tmp = Rivet::AnalysisLoader::getAllAnalysisNames();
        vector<string> knownAnalyses;
        for (set<string>::const_iterator i = tmp.begin(); i != tmp.end(); ++i) {
          knownAnalyses.push_back(*i);
          //knownAnalyses.push_back("~" + *i);
        }
        ValuesConstraint<string> anaNameConstraint(knownAnalyses);
        MultiArg<string> analysesArg("a", "analysis", "A Rivet analysis to be run. Prefix name with a ~ to disable instead", false, &anaNameConstraint, cmd);
        SwitchArg analysesAllArg("A", "all_analyses", "Run all Rivet analyses (unless any are disabled)", cmd, false);

        // Add histogram args
        ValueArg<string> histoNameArg("H", "histoname", "Base name of histogram file, with no extension ('Rivet' by default)", false, "Rivet", "name", cmd);
        vector<string> histformats = Rivet::getKnownHistoFormatNames();
        ValuesConstraint<string> histoTypeConstraint(histformats);
        ValueArg<string> histoTypeArg("", "histotype", "Histogram output format (default is AIDA XML)", false, "AIDA", &histoTypeConstraint, cmd);

        // Add logging args
        const string mesg = "Set log level in 'name=level' format. The levels are TRACE, DEBUG, INFO and WARN";
        MultiArg<string> logsArg("l", "loglevel", mesg, false, "logname=loglevel", cmd);
        SwitchArg disableLogColorArg("", "nocolor", "Disable shell color escapes (should be automatic for non-TTY stdout destination)", cmd, false);

        // Add misc args
        ValueArg<size_t> numEventsArg("n", "numevents", "Number of events to generate (10 by default)", false, 10, "num", cmd);
        ValueArg<string> hepmcOutFileArg("o", "outeventfile", "File to write HepMC events to (disabled by default)", false, "RivetGun.hepmc", "filename", cmd);
        SwitchArg disableRivetArg("R", "norivet", "Disable running of Rivet", cmd, false);

        /////////////////////////////////////////////////////////


        // Parse the command line args and make a Configuration
        cmd.parse(argc, argv);


        /////////////////////////////////////////////////////////


        // Have to read params first in case run params are included
        // Handle params from files/cin first...
        foreach (const string& pf, paramsFilesArg.getValue()) {
          if (pf != "-") {
            // cout << "Using param file " << pf << endl;
            std::ifstream in;
            in.open(pf.c_str());
            if (in.fail() && pf.find("/") == string::npos) {
              // Try finding the file in the standard installed RivetGun share directory
              const string try2 = getRivetgunDataPath() + "/" + pf;
              in.clear();
              in.open(try2.c_str());
              if (in.fail()) {
                throw Rivet::Error("Couldn't read from param file " + pf + " or " + try2);
              }
            }
            handleConfigStream(in, config.params);
            in.close();
          } else {
            handleConfigStream(cin, config.params);
          }
        }
        // ...then overload the file params with specific command line params
        foreach (const string& p, paramsArg.getValue()) {
          pair<string,string> key_val = parseParamString(p);
          if (!key_val.first.empty() && !key_val.second.empty()) {
            config.params.push_back(Param(key_val.first, key_val.second));
            //config.params[key_val.first] = key_val.second;
          } else {
            if (key_val.first.empty() && key_val.second.empty()) continue;
            throw Rivet::Error("Param setting key or value is empty, somehow: '" + p + "'");
          }
        }

        // Handle parameters which should be intercepted by RivetGun
        ParamMap::iterator p = config.params.begin();
        while (p != config.params.end()) {
          if (p->first.find("RG:") != string::npos) {
            //cout << "Found RG param: " << p->first << " = " << p->second << endl;
            stringstream pvalue(p->second);
            if (p->first == "RG:Mom1") {
              pvalue >> config.mom1;
            }
            else if (p->first == "RG:Mom2") {
              pvalue >> config.mom2;
            }
            else if (p->first == "RG:Beam1") {
              config.beam1 = Rivet::getParticleNameEnum(p->second);
            }
            else if (p->first == "RG:Beam2") {
              config.beam2 = Rivet::getParticleNameEnum(p->second);
            }
            else if (p->first == "RG:Seed") {
              pvalue >> config.rngSeed;
            }
            else {
              throw Rivet::Error("Unknown 'RG' parameter name: " + p->first);
            }

            // Remove this param so it doesn't get passed to the generator
            config.params.erase(p);
          } else {
            ++p;
          }
        }


        // Use generator args
        if (genNameArg.isSet()) {
          config.generatorName = genNameArg.getValue();
        } else if (genFileArg.isSet()) {
          /// @todo Read HepML file
          throw Rivet::Error("HepML file reading is not yet supported. Sorry.");
        }
        if (rngSeedArg.isSet()) config.rngSeed = rngSeedArg.getValue();
        if (eventFileArg.isSet()) {
          config.hepmcInFile = eventFileArg.getValue();
          config.readHepMC = true;
        }

        // Use initial state args
        try {
          if (beam1Arg.isSet()) config.beam1 = getParticleNameEnum(beam1Arg.getValue());
          if (beam2Arg.isSet()) config.beam2 = getParticleNameEnum(beam2Arg.getValue());
        } catch (exception& e) {
          throw Rivet::Error("Unknown beam particle: " + beam1Arg.getValue() + 
                              " or " + beam2Arg.getValue());
        }
        if (mom1Arg.isSet()) config.mom1 = mom1Arg.getValue();
        if (mom2Arg.isSet()) config.mom2 = mom2Arg.getValue();


        // Use analysis args
        // First handle the "enable all analyses" option...
        if (analysesAllArg.getValue()) {
          config.analyses = AnalysisLoader::getAllAnalysisNames();
        }
        // ...then handle individuals, including negations
        foreach (const string& a, analysesArg.getValue()) {
          string A = toUpper(a);
          try {
            // Check for analysis disabling with ~ prefix
            if (A.rfind("~", 0) == string::npos) {
              config.analyses.insert(A);
            } else {
              string Aneg = A.substr(1, A.size()-1);
              config.analyses.erase(Aneg);
            }
          } catch (std::exception& e) {
            throw Rivet::Error("Invalid analysis choice: " + a);
          }
        }


        // Use histo args
        // Get histogram filename
        config.histoName = histoNameArg.getValue();
        // Get histogram format
        if (histoTypeArg.getValue() == "AIDA") {
          config.histoFormat = AIDAML;
        } else if (histoTypeArg.getValue() == "FLAT") {
          config.histoFormat = FLAT;
        } else if (histoTypeArg.getValue() == "ROOT") {
          config.histoFormat = ROOT;
        }


        // Use logging args
        foreach (const string& l, logsArg.getValue()) {
          size_t breakpos = l.find("=");
          if (breakpos != string::npos) {
            string key = l.substr(0, breakpos);
            string value = l.substr(breakpos + 1, l.size() - breakpos - 1);
            const int level = Log::getLevelFromName(value);
            config.logLevels[key] = level;
          } else {
            throw Rivet::Error("Invalid log setting format: " + l);
          }
        }
        config.useLogColors = ! disableLogColorArg.getValue();


        // Misc. flags
        config.numEvents = numEventsArg.getValue();
        config.runRivet = ! disableRivetArg.getValue();
        if (hepmcOutFileArg.isSet()) {
          config.writeHepMC = true;
          config.hepmcOutFile = hepmcOutFileArg.getValue();
        }
      } catch (const TCLAP::ArgException& e) { 
        stringstream ss;
        ss << "Command line error: " << e.error() << " for arg " << e.argId(); 
        throw Rivet::Error(ss.str());
      }
      
      return config;
    }


  }
}


//////////////////////////////////////////////


namespace {
  using namespace std;


  pair<string,string> parseParamString(const string& paramstring) {
    pair<string,string> rtn;
    string line = paramstring;

    // Replace = signs and tabs with spaces
    while (line.find("=") != string::npos) {
      size_t eqpos = line.find("=");
      line.replace(eqpos, 1, " ");
    }
    while (line.find("\t") != string::npos) {
      size_t eqpos = line.find("\t");
      line.replace(eqpos, 1, " ");
    }
    /// @todo Strip spaces from ends with Boost string lib

    // If effectively empty, then end
    if (line.find_first_not_of(' ') == string::npos) {
      return rtn;
    }

    // If wrong format, say so with an Error
    if (line.find(' ') == string::npos) {
      throw Rivet::Error("Param setting string '" + paramstring + "' is not in a valid format.");
    }

    // Get name and value and add to the param stack by finding the
    // punctuation tokens at either end of the string.
    const size_t first = line.find_first_not_of(' ');
    const size_t last = line.find_last_not_of(' ');
    line = line.substr(first, last-first+1);
    //cout << "Stripped line :" << line << ":" << endl;
    const size_t firstgap = line.find(' ');
    const size_t lastgap = line.rfind(' ');
    const string pname = line.substr(first, firstgap-first);
    const string pval = line.substr(lastgap+1, last-lastgap);
    //cout << ":" << pname << "=" << pval << ":" << endl;
    
    // Check if the remaining bit of string isn't emptiness
    const string remains = line.substr(firstgap, lastgap-firstgap+1);
    if (remains.find_first_not_of(' ') != string::npos) {
      cerr << remains << " - error at char #" << remains.find_first_not_of(' ') << endl;
      throw Rivet::Error("Param setting string '" + paramstring + "' is not in a valid format.");
    }

    // If all is okay, set the parameter pair
    rtn.first = pname;
    rtn.second = pval;
    return rtn;
  }
  
  
  void handleConfigStream(istream& in, Rivet::ParamMap& pmap) {
    while(in) {
      string line;
      getline(in, line);
      const string origline = line;
      //cout << "Input param line: " << line << endl;

      // Chop line after a # character
      size_t commentPosn = line.find_first_of("#");
      if (commentPosn != string::npos) {
        line = line.substr(0, commentPosn);
      }

      // Parse param string
      pair<string,string> key_val = parseParamString(line);
      if (!key_val.first.empty() && !key_val.second.empty()) {
        pmap.push_back(Rivet::Param(key_val.first, key_val.second));
        //pmap[key_val.first] = key_val.second;
        //if(key_val.first.find("RG:") == string::npos) nmap.push_back(key_val.first);
      } else {
        continue;
        //throw Rivet::Error("Param setting key or value is empty, somehow: '" + line + "'");
      }

    }


  }
}
