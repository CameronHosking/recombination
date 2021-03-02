#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <unordered_map>
#include <vector>
#include <algorithm>

bool startMatches(const std::string &a,const std::string &b)
{
    for(uint i = 0 ;i < std::min(a.size(),b.size());++i)
    {
        if (a[i]!=b[i])
        {
            return false;
        }
    }
    return true;
}
enum MutationType{Valid, Ignored, Fatal};
std::unordered_map<int,int> ignoredSiteCounts;
MutationType validMutation(const std::string &mutation)
{
    //check that the sequence contains no ambiguity characters
    uint lastTab = mutation.find_last_of('\t');

    std::string site = mutation.substr(lastTab+1);
    site = site.substr(0,site.find(':'));
    int siteNumber = std::stoi(site);
    //the mutation is near the start or end of the sequence
    if(siteNumber < 50||siteNumber>29750)
    {
        ignoredSiteCounts[siteNumber]++;
        return Ignored;
    }
    //the first non Nucleotide character is within the ref/alt columns
    //
    if(mutation.find_first_not_of("acgtACGT-\t")<lastTab)
    {
        return Fatal;
    }

    return Valid;
    /*
    //if at the start or the end
    if(site.substr(0,2)=="1:")
        return false;
    if(site.find(":29902")!=std::string::npos)
        return false;
    //check for dissallowed positions
    if(startMatches(site,"635:")) return false;
    if(startMatches(site,"2091:")) return false;
    if(startMatches(site,"2094:")) return false;
    if(startMatches(site,"3145:")) return false;
    if(startMatches(site,"3564:")) return false;
    if(startMatches(site,"4050:")) return false;
    if(startMatches(site,"5736:")) return false;
    if(startMatches(site,"6869:")) return false;
    if(startMatches(site,"8022:")) return false;
    if(startMatches(site,"8790:")) return false;
    if(startMatches(site,"10129:")) return false;
    if(startMatches(site,"11074:")) return false;
    if(startMatches(site,"11083:")) return false;
    if(startMatches(site,"11535:")) return false;
    if(startMatches(site,"13402:")) return false;
    if(startMatches(site,"13408:")) return false;
    if(startMatches(site,"13476:")) return false;
    if(startMatches(site,"13571:")) return false;
    if(startMatches(site,"14277:")) return false;
    if(startMatches(site,"15922:")) return false;
    if(startMatches(site,"16887:")) return false;
    if(startMatches(site,"19484:")) return false;
    if(startMatches(site,"21575:")) return false;
    if(startMatches(site,"22335:")) return false;
    if(startMatches(site,"24389:")) return false;
    if(startMatches(site,"24390:")) return false;
    if(startMatches(site,"24933:")) return false;
    if(startMatches(site,"26549:")) return false;
    if(startMatches(site,"29037:")) return false;
    if(startMatches(site,"29553:")) return false;
    */    
}

//counts the number of elements in common, the elements must be sorted
int elementsInCommon(const std::vector<int> &a, const std::vector<int> &b, const std::vector<int> &c)
{
    uint i = 0;
    uint j = 0;
    uint k = 0;
    int common = 0;
    while(i < a.size() && j < b.size() && k < c.size())
    {
        if(a[i]==b[j]&&a[i]==c[k])
        {
            common++;
            i++;
            j++;
            k++;
        }
        else if(a[i]>b[j])
        {
            if(c[k]>b[j])
                j++;
            else
                k++;         
        }
        else
        {
            if(c[k]>a[i])
                i++;
            else
                k++;  
        }
    }
    return common;
}
//map of mutation string to mutation ID
std::unordered_map<std::string,int> mutationIDs;

//vector of strains, this ordering represents their numeric key
std::vector<std::string> strains;

//vector of strains by strainID that contain each mutation by mutationID
std::vector<std::vector<int> > strainsContainingMutation;

//vector of the mutations each strain contains
std::vector<std::vector<int> > strainMutations;

std::unordered_map<std::string,std::string> uniqueStrainMutations;

void deleteLastStrainAddedIfStrainIs(const std::string& strain)
{
    if(!strains.empty()&&strains.back()==strain)
    {
        int strainID = strains.size()-1;
        //remove the strain from list of strains
        strains.pop_back();
        //remove this strains from all the mutations that have this strain
        for(auto m: strainMutations[strainID])
        {
            strainsContainingMutation[m].pop_back();
        }
        //remove this strain and its mutations from the list of strains and their mutations
        strainMutations.pop_back();
    }
};

int main(int argc, char* argv[])
{
     // Check the number of parameters
    if (argc < 4) {
        // Tell the user how to run the program
        std::cerr << "Usage: " << argv[0] << " INPUTFILE OUTPUTFILE DUPLICATESFILE" << std::endl;
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }


    std::ifstream infile(argv[1]);
    std::ofstream duplicatesFile(argv[3]);
    duplicatesFile << "Original Duplicate" << std::endl;

    std::string line;



    std::cout << "reading in mutations from " << argv[1] << std::endl;
    std::cout << "writing duplicates to " << argv[3] << std::endl;
    //this strain contained ambiguity codes
    bool invalidStrain = false;
    std::string prevStrain = "";
    std::string allMutationsString = "";
    while (std::getline(infile, line))
    {
        int firstTab = line.find('\t');
        std::string strain = line.substr(0,firstTab);
        if(strain!=prevStrain)
        {
            if(!invalidStrain)
            {
                //if this is a new strain chack the last strain to make sure its mutation set was unique
                if(uniqueStrainMutations.count(allMutationsString)==0)
                {
                    uniqueStrainMutations[allMutationsString]=prevStrain;
                }
                else
                {
                    deleteLastStrainAddedIfStrainIs(prevStrain);
                    duplicatesFile << prevStrain << "\t" << uniqueStrainMutations[allMutationsString] <<std::endl;
                }
            }
            //reset mutation string
            allMutationsString = "";
            invalidStrain = false;
            prevStrain = strain;
        }
        else
        {
            //if this strain has had invalid mutations ignore it
            if(invalidStrain)
            {
                continue;
            }
        }
        
        std::string mutation = line.substr(firstTab+1);
        //check that mutation is valid
        MutationType mutationType = validMutation(mutation);
        if(mutationType==Fatal)
        {
            invalidStrain = true;
            //remove all trace of this strain
            deleteLastStrainAddedIfStrainIs(strain);
            continue;
        }

        if(mutationType==Ignored)
            continue;

        //valid muation is added to this strains all mutations string
        allMutationsString += mutation;
        //if it is a new strain we add it to the list of strains
        if(strains.empty()||strains.back()!=strain)
        {
            strains.push_back(strain);
            strainMutations.push_back(std::vector<int>());
        }
        int strainID = strains.size()-1;
        
        //if its a new mutation add it to the mutation set
        if(mutationIDs.count(mutation)==0)
        {
            mutationIDs[mutation] = mutationIDs.size();
            strainsContainingMutation.push_back(std::vector<int>());
        }
        int mutationID = mutationIDs[mutation];
        
        //add this strain to the list of strains that have this mutation
        strainsContainingMutation[mutationID].push_back(strainID);
        //and add this mutation to this strains list of mutations
        strainMutations[strainID].push_back(mutationID);
    }
    duplicatesFile.close();
    std::cout << "read in " << mutationIDs.size() << " unique mutations in " << strains.size() << " strains" << std::endl;
    
    // std::vector<std::pair<int,int>> sitesIgnored = std::vector<std::pair<int,int>>(ignoredSiteCounts.begin(),ignoredSiteCounts.end());
    // std::cout << sitesIgnored.size() <<std::endl;
    // std::sort(sitesIgnored.begin(),sitesIgnored.end());
    // for(auto s: sitesIgnored)
    // {
    //     std::cout << s.first << "\t" << s.second << std::endl;
    // }
    //return 0;

    //sort the mutations (this will be useful later when comparing mutations between strains)
    for(uint i = 0; i < strains.size();i++)
    {
        std::sort(strainMutations[i].begin(),strainMutations[i].end());
    }

    std::cout << "creating common mutatations matrix" << std::endl;
    //create commonMutations matrix, initially zeroes
    std::vector<std::vector<unsigned char> > commonMutations;
    for(uint i = 0; i < strains.size(); ++i)
    {
        commonMutations.push_back(std::vector<unsigned char>(i+1,0));
    }

    auto numberOfCommonMutations = [&](int strainA, int strainB) -> unsigned char&
    {
        //swap strainA and strainB if A is smaller than B, that way
        //we only have to fill out and access half of the AB matrix
        if(strainA<strainB)
        {
            std::swap(strainA,strainB);
        }
        return commonMutations[strainA][strainB];
    };

    //loop through each mutation, adding one to the count of mutations in common to each strain pair
    for(uint i = 0; i < strainsContainingMutation.size(); i++)
    {
        for(uint s1 = 0; s1 < strainsContainingMutation[i].size(); s1++)
        {                
            int s1ID = strainsContainingMutation[i][s1];
            for(uint s2 = 0; s2 <= s1; s2++)
            {
                int s2ID = strainsContainingMutation[i][s2];
                numberOfCommonMutations(s1ID,s2ID)++;
            }
        }
    }
    std::cout << "finished creating common mutatations matrix" << std::endl;

    int maxInCommon = 0;
    for(uint i = 0; i < strains.size(); ++i)
    {
        for(uint j = 0; j < i; ++j)
        {
            int inCommon = commonMutations[i][j];
            if(inCommon>maxInCommon)
                maxInCommon = inCommon;
        }
    }
    //data structure to store strains ordered by similarity
    std::vector<std::vector<int> > strainsBySimilarity = std::vector<std::vector<int> >(maxInCommon+1,std::vector<int>());

    std::ofstream outfile(argv[2]);


    // struct int30{
    //     int nums[30];
    // };
    // int30 distancesToNearestParent[60] = {};
    // int10 differentParents[60] = {};
    for(uint i = 0; i < strains.size(); ++i)
    {
        if(i<38904) continue;
        if(i%100==0)
        {
            std::cout << "Tested " << i << "/"<< strains.size() <<" sequences" << std:: endl;
        }
        //check to see if the strain has a direct parent
        int iMutations = commonMutations[i][i];
        int bestDirectParent = i;
        int maxAccountedFor = 0;
        bool directParentFound = false;
        //int duplicateParents = 0;
        for(uint j = 0; j < strains.size(); j++)
        {
            int jMutations = commonMutations[j][j];
            int numijCommonMutations = numberOfCommonMutations(i,j);

            //i is completely a subset of j mutations and is hence an ancestor of j (or identical)
            if(numijCommonMutations==iMutations)
                continue;

            //not all of js mutations are in common with i so not a direct parent
            if(jMutations != numijCommonMutations)
                continue;
            
            //if j has the most mutations in common seen so far but still fewer than i
            //it is the new best direct parent
            if(numijCommonMutations > maxAccountedFor)
            {
                //duplicateParents = 0;
                bestDirectParent = j;
                maxAccountedFor = numijCommonMutations;
                //direct parent within distance 2 found
                if(numijCommonMutations==iMutations-1||numijCommonMutations==iMutations-2)
                {
                    directParentFound = true;
                    break;
                }
            }
            // else if(numijCommonMutations == maxAccountedFor)
            // {
            //     duplicateParents++;
            // }           
        }
        // int dist = iMutations - maxAccountedFor;
        // if(iMutations < 30 && dist < 20)
        //     distancesToNearestParent[iMutations].nums[dist]++;
        // continue;
        // if(duplicateParents>9)
        // {
        //     duplicateParents = 9;
        // }
        // differentParents[iMutations - maxAccountedFor].nums[duplicateParents]++;
        // continue;

        //don't need to look for recombination events if we've already got a direct parent
        if(directParentFound){
            continue;
        }

        //clear the strains by similarity data structure to fill it with the new strain info
        for(uint k = 0; k < strainsBySimilarity.size();k++)
        {
            strainsBySimilarity[k].clear();
        }
        //fill out the strains by similarity data structure, now strains are sorted by similarity.
        int maxSimilarity = 0;
        for(uint j = 0; j < strains.size(); j++)
        {
            int numijCommonMutations = numberOfCommonMutations(i,j);

            //if the strains is sufficiently similar but not likely to be a descendant store it
            //it could potentially be a parent
            if(numijCommonMutations>1 && numijCommonMutations<iMutations-1)
            {
                maxSimilarity = std::max(maxSimilarity,numijCommonMutations);
                strainsBySimilarity[numijCommonMutations].push_back(j);
            }
        }

        //check for recombination events
        int minExtraMutations = 0;
        int bestP1 = i;
        int bestP2 = i;
        //int64_t comparisons = 0;
        //similarity must be greater than or equal to 3 for it to be possible to have 3 unique
        for(int p1Similarity = maxSimilarity;p1Similarity >=3; p1Similarity--)
        {
            for(auto p1It = strainsBySimilarity[p1Similarity].begin();p1It!=strainsBySimilarity[p1Similarity].end();p1It++)
            {
                int p1 = *p1It;
                int mutationsInP1 = commonMutations[p1][p1];

                //keep checking against parents with lower and lower similarity until it's impossible that they could be an improvement
                for(int p2Similarity = p1Similarity; p2Similarity+p1Similarity >= maxAccountedFor&&p2Similarity>=3; p2Similarity--)
                {
                    //loop through all strains with similarity p2, if p1 and p2 are the same we only want to loop on sequences after p1 to avoid duplicating work
                    for(auto p2It = (p1Similarity==p2Similarity) ? p1It+1 : strainsBySimilarity[p2Similarity].begin(); p2It!=strainsBySimilarity[p2Similarity].end();p2It++)
                    {
                        //comparisons++;
                        int p2=*p2It;
                        //mutations accounted for, assumes parent common mutations are in child
                        int accountedFor = p1Similarity + p2Similarity - numberOfCommonMutations(p1,p2);

                        //we have a better parent / parents
                        if(accountedFor < maxAccountedFor)
                        continue;

                        int mutationsInP2 = commonMutations[p2][p2];
                        int extraMutations = mutationsInP1 - p1Similarity + mutationsInP2 - p2Similarity;
                        
                        //pick parents with the minimum number of extraneous mutations
                        if(accountedFor==maxAccountedFor&&extraMutations>minExtraMutations)
                            continue;

                        //must have at least 3 unique mutations from each parent, this is to avoid the case where a recombination occurs but we have only sequenced the children
                        //of the recombinant strain and each selects their sibling as their recombination parent
                        int p1p2SharedMutations = numberOfCommonMutations(p1,p2);
                        if(p1Similarity-p1p2SharedMutations < 3 || p2Similarity-p1p2SharedMutations < 3)
                            continue;
                        //make sure that the mutations in common between parents are actually present in child
                        if(elementsInCommon(strainMutations[i],strainMutations[p1],strainMutations[p2])==numberOfCommonMutations(p1,p2))
                        {
                            //if we get this far this is our new best guess for parents
                            maxAccountedFor = accountedFor;
                            minExtraMutations = extraMutations;
                            bestP1 = p1;
                            bestP2 = p2;
                        }
                    }
                } 
            }
        }
        //std::cout << comparisons << std::endl;
        //at least 3 better than a direct parent
        if(maxAccountedFor-3>=numberOfCommonMutations(i,bestDirectParent))
        {
            int commonWithP1 = numberOfCommonMutations(i,bestP1);
            int commonWithP2 = numberOfCommonMutations(i,bestP2);
            //difference between recombination event and direct parent in number of mutations explained
            outfile << maxAccountedFor - numberOfCommonMutations(i,bestDirectParent) << '\t';
            // //minimum unique
            // int minimumUnique = maxAccountedFor - std::max(commonWithP1,commonWithP2);
            // if(minimumUnique < 2)
            //     continue;
            // outfile << minimumUnique << '\t';
            
            //direct parent number of mutations explained
            outfile << (int) (commonMutations[i][i] - numberOfCommonMutations(i,bestDirectParent)) << '\t';
            //distribution of mutations: mutations unique to p1 / number of shared mutations / mutations unique to p2
            int shared = numberOfCommonMutations(bestP1,bestP2);
            outfile << commonWithP1 - shared << '/' << shared << '/' << commonWithP2 - shared << '\t';
            //mutations transfered in event compared to total mutations in strain
            outfile << maxAccountedFor << '/' << (int) commonMutations[i][i] <<'\t';
            outfile << commonWithP1 << '/' << (int) commonMutations[bestP1][bestP1] << '\t';
            outfile << commonWithP2 << '/' << (int) commonMutations[bestP2][bestP2] << '\t';
            //actual strain names
            outfile << strains[i] << '\t';
            outfile << strains[bestP1] << '\t';
            outfile << strains[bestP2] << std::endl;
        }
    }
    // for(uint j = 1; j < 30; j++)
    // {
    //     std::cout << '\t' << j;
    // }
    // std::cout << std::endl;
    // for(uint i = 1; i < 30; i++)
    // {
    //     std::cout << i;
    //     for(uint j = 1; j < 20; j++)
    //     {
    //         std::cout << '\t' << distancesToNearestParent[i].nums[j];
    //     } 
    //     std::cout << std::endl;
    // }
    return 0;
}
