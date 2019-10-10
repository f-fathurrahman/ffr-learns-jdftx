#include <iostream>
#include <core/string.h>
#include <limits.h>
#include <algorithm>

#include "jdftx/commands/command.h"
#include "jdftx/electronic/Everything.h"
#include "jdftx/commands/parser.h"
#include "jdftx/electronic/ColumnBundle.h"

using std::map;
using std::pair;
using std::make_pair;
using std::set;
using std::cout;
using std::endl;

//! Process the statically initialized command map, add dependency info
//! and check for inconsistencies (static ones i.e. bugs in the code, not input file)
struct ProcessedCommandMap : public map<string, pair<int,Command*> >
{
    int nPasses; //!< maximum number of passes required based on dependency tree

    //! Calls getCommandMap(), checks and encodes dependency information
    //! commands paired with integer 0 have no dependencies,
    //! those with 1 depend only on others with 0 etc.
    ProcessedCommandMap()
    {
        cout << "Pass here 26" << endl;
        const map<string,Command*>& origMap = getCommandMap();
        //Copy into *this with all the integers set to INT_MAX
        for(map<string,Command*>::const_iterator i=origMap.begin(); i!=origMap.end(); i++)
            (*this)[i->first] = make_pair(INT_MAX, i->second);
        //Make multiple passes till all the commands have been handled:
        for(nPasses=0; ; nPasses++)
        {   
            cout << "nPasses = " << nPasses << endl;
            int nProcessed = 0; //number processed in this pass
            int nRemaining = 0; //number not yet processed
            for(iterator i=begin(); i!=end(); i++)
                if(i->second.first>nPasses) //not yet processed
                { 
                    Command& ci = *(i->second.second);
                    bool ready = true;
                    cout << "Pass here 45" << endl;
                    //Loop through all its required commands:
                    for(set<string>::iterator r=ci.requires.begin(); r!=ci.requires.end(); r++)
                    {
                        cout << "Pass here 49" << endl;
                        iterator j = find(*r);
                        if(j==end())
                            die("Inconsistency in command map: '%s' in '%s's requires list is not a valid command.\n",
                                r->c_str(), ci.name.c_str());
                        if(j->second.first>=nPasses) { ready=false; break; } //unsatisfied dependency
                    }
                    cout << "Pass here 56" << endl;
                    if(ready)
                    {   i->second.first = nPasses;
                        nProcessed++;
                    }
                    else nRemaining++;
                }
            if(!nRemaining) break; //all commands have been processed
            if(!nProcessed)
            {   logPrintf("There is a cyclic dependency somewhere amongst the following commands:\n");
                for(iterator i=begin(); i!=end(); i++)
                    if(i->second.first>nPasses)
                        logPrintf("\t%s\n", i->second.second->name.c_str());
                die("Cyclic dependency in command map.\n");
            }
        }
        //Make forbids mutual:
        typedef map<string,Command*>::const_iterator Citer;
        for(Citer i=origMap.begin(); i!=origMap.end(); i++)
        {   Command& ci = *(i->second);
            for(set<string>::iterator f=ci.forbids.begin(); f!=ci.forbids.end(); f++)
            {   Citer j = origMap.find(*f);
                if(j==origMap.end())
                    die("Inconsistency in command map: '%s' in '%s's forbids list is not a valid command.\n",
                        f->c_str(), ci.name.c_str());
                Command& cj = *(j->second);
                if(!cj.forbids.count(ci.name) && cj.section==ci.section) //don't report cross section errors as they can't be fixed
                {   logPrintf("Command map WARNING: '%s' forbids '%s', but not vice versa (making this so)\n",
                        ci.name.c_str(), cj.name.c_str());
                    cj.forbids.insert(ci.name);
                }
            }
        }
        //Check for inconsistency in the forbids:
        for(Citer i=origMap.begin(); i!=origMap.end(); i++)
        {   Command& ci = *(i->second);
            set<Command*> deps, forbids;
            getDependencies(ci, origMap, deps, forbids);
            //Find the intersection of the dependencies and the forbidden commands:
            std::vector<Command*> common(std::max(deps.size(),forbids.size()));
            common.erase(
                set_intersection(deps.begin(),deps.end(), forbids.begin(),forbids.end(), common.begin()),
                common.end() );
            if(common.size())
            {   logPrintf("Command '%s' (in)directly requires as well as forbids the following commands:\n", ci.name.c_str());
                for(unsigned j=0; j<common.size(); j++)
                    logPrintf("\t%s\n", common[j]->name.c_str());
                die("Forbid/require inconsistency in command map.\n");
            }
        }
    }

private:
    //Add all the dependencies of command c in dep,
    //and all the commands forbidden by c and its dependencies to forbids
    void getDependencies(const Command& c, const map<string,Command*>& cmdMap, set<Command*>& dep, set<Command*>& forbids)
    {   for(set<string>::iterator r=c.requires.begin(); r!=c.requires.end(); r++)
        {   Command* cr = cmdMap.find(*r)->second;
            dep.insert(cr);
            for(set<string>::iterator f=cr->forbids.begin(); f!=cr->forbids.end(); f++)
                forbids.insert(cmdMap.find(*f)->second);
            getDependencies(*cr, cmdMap, dep, forbids);
        }
    }
};


//! Call Command::process with error handling, count updating etc:
void safeProcess(Command& c, string params, Everything& everything,
    map<string,int>& encountered, std::vector< pair<Command*,string> >& errors);

int main(int argc, char** argv)
{

  //std::vector< pair<string,string> > input, Everything& everything, bool printDefaults

 //Parse command line, initialize system and logs:
  Everything everything; //the parent data structure for, well, everything
  
  InitParams ip("Performs Joint Density Functional Theory calculations.", &everything);
  initSystemCmdline(argc, argv, ip);
  
  //Parse input file and setup
  ElecVars& eVars = everything.eVars;

  cout << "Filename = " << ip.inputFilename << endl;

  auto input = readInputFile(ip.inputFilename);
  auto printDefaults = ip.printDefaults;

  ProcessedCommandMap cmap;

  set<string> unknown; //unknown command names
  map<string,int> encountered; //command names, and the number of times they were encountered
  std::vector< pair<Command*,string> > errors; //command classes that encountered errors, with error messages
  std::vector<std::shared_ptr<SpeciesInfo>> species; //temporary copy of ion-info species

//First check for, take note of, and remove unknown commands
    for(unsigned i=0; i<input.size();)
    { 
        cout << "i = " << i << endl;
        if(cmap.find(input[i].first)==cmap.end())
        {   unknown.insert(input[i].first);
            input.erase(input.begin()+i);
        }
        else i++;
    }
    
    //All command names in input are now known to be in cmap (so can safely use cmap[] instead of cmap.find())
    for(int pass=0; pass<=cmap.nPasses; pass++)
    {   
        cout << "Pass = " << pass << endl;
        //Run through all commands marked for this pass:
        bool encounteredIon = false;
        for(unsigned i=0; i<input.size(); i++)
        {   pair<int,Command*>& icPair = cmap[input[i].first];
            if(icPair.first==pass)
            {   Command& c = *(icPair.second);
                if(c.name == "ion") encounteredIon = true;
                //Check for empty parameter list:
                string trimmedParams = input[i].second;
                trim(trimmedParams);
                if(!trimmedParams.length() && c.emptyParamError.length())
                    errors.push_back(make_pair(&c,
                        " cannot be issued without any parameters:\n"
                        + c.emptyParamError
                        + (c.hasDefault ? "\nOmit this command altogether for the default behaviour.\n" : "")));
                //Process command:
                string params = input[i].second + " "; //add space at end to prevent eof on last parameter
                safeProcess(c, params, everything, encountered, errors);
            }
        }
        //Run defaults where necessary:
        for(ProcessedCommandMap::iterator i=cmap.begin(); i!=cmap.end(); i++)
            if(i->second.first==pass)
            {   Command& c = *(i->second.second);
                if(c.hasDefault && !encountered[c.name])
                    safeProcess(c, "", everything, encountered, errors);
            }
        //Remove unused pseudopotentials as soon as ion command has been processed:
        if(encounteredIon)
        {   species = everything.iInfo.species; //backup of species that retains unused ones
            for(auto iter=everything.iInfo.species.begin(); iter!=everything.iInfo.species.end();)
                if(not (*iter)->atpos.size())
                    iter = everything.iInfo.species.erase(iter);
                else iter++;
        }
    }
    //Quit on errors:
    size_t errTot = unknown.size() + errors.size();
    if(errTot)
    {   logPrintf("\n\nEncountered %lu errors while parsing input:\n\n", errTot);
        for(set<string>::iterator i=unknown.begin(); i!=unknown.end(); i++)
            logPrintf("'%s' is not a valid command.\n", i->c_str());
        logPrintf("\n");
        for(unsigned i=0; i<errors.size(); i++)
            logPrintf("Command '%s' %s\n\n", errors[i].first->name.c_str(), errors[i].second.c_str());
        die("\n\nInput parsing failed with %lu errors (run with -t for command syntax)\n\n", errTot);
    }
    //Print status:
    std::swap(species, everything.iInfo.species); //present original list of species for printStatus
    std::map<string,int> explicitlyEncountered; //List of commands explicitly in input file (with corresponding multiplicities)
    if(!printDefaults) for(auto cmd: input) explicitlyEncountered[cmd.first]++;
    logPrintf("\n\nInput parsed successfully to the following command list (%sincluding defaults):\n\n", printDefaults ? "" : "not ");
    for(auto i: (printDefaults ? encountered : explicitlyEncountered))
    {   for(int iRep=0; iRep<i.second; iRep++) //handle repetitions
        {   Command& c = *(cmap[i.first].second);
            logPrintf("%s ", c.name.c_str());
            c.printStatus(everything, iRep);
            logPrintf("\n");
        }
    }
    logPrintf("\n\n");
    std::swap(species, everything.iInfo.species); //restore list of species with unused ones removed

//  cout << "species = " << species << endl;
  cout << "Pass here ..." << endl;
  return 0;
}
