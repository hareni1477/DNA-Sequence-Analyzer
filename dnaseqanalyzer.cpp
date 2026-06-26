#include <bits/stdc++.h>
#include <cctype>
using namespace std;
vector<string> orfs; // Stores all Open Reading Frames (ORFs) found in the DNA sequence.
// Searches the DNA sequence for all occurrences of a user-specified motif
void motifSearch(const string &dna,string motif)
{
    int n = motif.length();
    for(char &ch : motif)
{
    ch = toupper(ch);
}

    bool found = false;
    for(int i = 0; i <= dna.length() - n; i++){
        if(dna.substr(i, n) == motif)
    {
        cout << "Motif found at position "
             << i + 1 << endl;

        found = true;
    }
    }
    if(!found)
{
    cout << "Motif not found in the DNA sequence." << endl;
}
}   
// Finds all valid Open Reading Frames (ORFs).
// An ORF starts with ATG and ends with the first in-frame stop codon (TAA, TAG or TGA).
void findORFs(const string &dna)
{
    
    orfs.clear();    
    bool foundORF = false;

    for (int i = 0; i <= dna.length() - 3; i++)
    {
        string codon = dna.substr(i, 3);

        if (codon == "ATG")
        {
            bool stopFound = false;

            // Search only in the same reading frame
            for (int j = i + 3; j <= dna.length() - 3; j += 3)
            {
                string stopCodon = dna.substr(j, 3);

                if (stopCodon == "TAA" ||
                    stopCodon == "TAG" ||
                    stopCodon == "TGA")
                {
                    foundORF = true;
                    stopFound = true;

                    cout << "\nORF Found!" << endl;
                    cout << "Start Position : " << i + 1 << endl;
                    cout << "Stop Position  : " << j + 3 << endl;

                    string orf = dna.substr(i, j + 3 - i);

                    orfs.push_back(orf); // Store the ORF in the global vector

                    cout << "ORF Sequence   : "
                         << orf << endl;

                    break;
                }
            }

            if (!stopFound)
            {
                cout << "\nStart codon found at position "
                     << i + 1
                     << " but no valid stop codon found."
                     << endl;
            }
        }
    }

    if (!foundORF)
    {
        cout << "\nNo ORFs found in the DNA sequence."
             << endl;
    }
}
// Standard genetic code used to translate DNA codons into amino acids.
unordered_map<string, char> geneticCode =
{
    // Phenylalanine
    {"TTT",'F'}, {"TTC",'F'},

    // Leucine
    {"TTA",'L'}, {"TTG",'L'},
    {"CTT",'L'}, {"CTC",'L'},
    {"CTA",'L'}, {"CTG",'L'},

    // Isoleucine
    {"ATT",'I'}, {"ATC",'I'}, {"ATA",'I'},

    // Methionine (Start)
    {"ATG",'M'},

    // Valine
    {"GTT",'V'}, {"GTC",'V'},
    {"GTA",'V'}, {"GTG",'V'},

    // Serine
    {"TCT",'S'}, {"TCC",'S'},
    {"TCA",'S'}, {"TCG",'S'},
    {"AGT",'S'}, {"AGC",'S'},

    // Proline
    {"CCT",'P'}, {"CCC",'P'},
    {"CCA",'P'}, {"CCG",'P'},

    // Threonine
    {"ACT",'T'}, {"ACC",'T'},
    {"ACA",'T'}, {"ACG",'T'},

    // Alanine
    {"GCT",'A'}, {"GCC",'A'},
    {"GCA",'A'}, {"GCG",'A'},

    // Tyrosine
    {"TAT",'Y'}, {"TAC",'Y'},

    // Histidine
    {"CAT",'H'}, {"CAC",'H'},

    // Glutamine
    {"CAA",'Q'}, {"CAG",'Q'},

    // Asparagine
    {"AAT",'N'}, {"AAC",'N'},

    // Lysine
    {"AAA",'K'}, {"AAG",'K'},

    // Aspartic Acid
    {"GAT",'D'}, {"GAC",'D'},

    // Glutamic Acid
    {"GAA",'E'}, {"GAG",'E'},

    // Cysteine
    {"TGT",'C'}, {"TGC",'C'},

    // Tryptophan
    {"TGG",'W'},

    // Arginine
    {"CGT",'R'}, {"CGC",'R'},
    {"CGA",'R'}, {"CGG",'R'},
    {"AGA",'R'}, {"AGG",'R'},

    // Glycine
    {"GGT",'G'}, {"GGC",'G'},
    {"GGA",'G'}, {"GGG",'G'},

    // Stop Codons
    {"TAA",'*'},
    {"TAG",'*'},
    {"TGA",'*'}
};
string currentProtein = "";
void translateORF(const vector<string> &orfs)
{
    
    if (orfs.empty())
    {
        cout << "No ORFs available for translation.Run findorf first." << endl;
        return;
    }
    else {
        cout << "Choose an ORF to translate ";
        for (int i = 0; i < orfs.size(); i++)
        {
            cout << "\n" << i + 1 << ". " << orfs[i];
        }
        cout << "\nEnter your choice: ";
        int choice;
        cin >> choice;
        
        string protein = "";
        if (choice < 1 || choice > orfs.size())
        {
            cout << "Invalid choice!" << endl;
            return;
        }
        string selectedORF = orfs[choice - 1];
    
        for (int i = 0; i <= selectedORF.length() - 3; i += 3)
        {
            string codon = selectedORF.substr(i, 3);
            if (geneticCode.find(codon) == geneticCode.end())
            {
                cout << "Unknown codon: " << codon << endl;
                return;
            }
            char aminoAcid = geneticCode[codon];
            if (aminoAcid == '*')
            break;
            protein += aminoAcid;
        }
        currentProtein = protein;
        cout << "Protein Sequence: " << protein << endl;
    }
}
unordered_map<char, double> aminoAcidWeights = {
        {'A', 89.09}, {'R', 174.20}, {'N', 132.12}, {'D', 133.10},
        {'C', 121.15}, {'E', 147.13}, {'Q', 146.15}, {'G', 75.07},
        {'H', 155.16}, {'I', 131.17}, {'L', 131.17}, {'K', 146.19},
        {'M', 149.21}, {'F', 165.19}, {'P', 115.13}, {'S', 105.09},
        {'T', 119.12}, {'W', 204.23}, {'Y', 181.19}, {'V', 117.15}
    };


void proteinstatisticsandweight(const string &currentProtein)
{
    if (currentProtein.empty())
    {
        cout << "No protein sequence available. Please translate an ORF first." << endl;
        return;
    }
    unordered_map<char, int> aminoAcidCount;

    
    double totalWeight = 0;
    for (char aa : currentProtein)
    {   
        auto it = aminoAcidWeights.find(aa);
        if(it != aminoAcidWeights.end())
        {
            aminoAcidCount[aa]++;
            totalWeight += it->second;
        }
        else
        {
            cout << "Unknown amino acid: " << aa << endl;
            return;
        }
        
    }

    cout << "Protein Sequence: " << currentProtein << endl;
    cout << "Length: " << currentProtein.length() << " amino acids" << endl;
    cout << "Molecular Weight: " << totalWeight << " Da" << endl;
    for(auto pair : aminoAcidCount)
    {
        cout << pair.first
         << " : "
         << pair.second
         << "\n";
    }
}
int main()
{
    string dna_sequence;

    cout << "Enter a DNA sequence: ";
    cin >> dna_sequence;

    int countA = 0;
    int countT = 0;
    int countC = 0;
    int countG = 0;
    string reverse_complement = "";
    // Validation + Uppercase Conversion + Counting
    for (char &ch : dna_sequence)
    {
        ch = toupper(ch);

        if (ch == 'A')
        {
            countA++;
            reverse_complement = 'T' + reverse_complement;
        }
        else if (ch == 'T')
        {
            countT++;
            reverse_complement = 'A' + reverse_complement;
        }
        else if (ch == 'C')
        {
            countC++;
            reverse_complement = 'G' + reverse_complement;
        }
        else if (ch == 'G')
        {
            countG++;
            reverse_complement = 'C' + reverse_complement;
        }
        else
        {
            cout << "Invalid character in DNA sequence: "
                 << ch << endl;
            return 1;
        }
    }

    int choice;

    while (true)
    {
        cout << "\n===== DNA SEQUENCE ANALYZER =====\n";

cout << "\n--- Sequence Analysis ---\n";
cout << "1. Show DNA Sequence\n";
cout << "2. Validate DNA Sequence\n";
cout << "3. Count Nucleotides\n";
cout << "4. Calculate GC Content\n";
cout << "5. Show Reverse Complement\n";

cout << "\n--- Gene Analysis ---\n";
cout << "6. Search Motif\n";
cout << "7. Find ORFs\n";

cout << "\n--- Protein Analysis ---\n";
cout << "8. Translate ORF\n";
cout << "9. Protein statistics and molecular weight\n";
cout << "\n10. Exit\n";

cout << "Enter choice: ";

        cin >> choice;

        switch (choice)
        {
        case 1:
        
            cout << "DNA Sequence: "
                 << dna_sequence << endl;
            break;
        


        case 2:
        
            cout << "DNA Sequence is Valid" << endl;
            break;
        
        case 3:
        
            cout << "A: " << countA << endl;
            cout << "T: " << countT << endl;
            cout << "C: " << countC << endl;
            cout << "G: " << countG << endl;
            break;
    

        case 4:
        {
            double gc_content =
                ((double)(countG + countC) * 100) /
                dna_sequence.length();

            cout << "GC Content: "
                 << gc_content
                 << "%" << endl;

            break;
        }

        case 5:
            cout << "Reverse Complement: "
                 << reverse_complement << endl;
            break;
        
        case 6:
        {
            string motif;
            cout << "Enter motif : ";
            cin >> motif;
            motifSearch(dna_sequence, motif);
            break;
        }
        case 7:
        {
            cout << "Finding ORFs..." << endl;
            findORFs(dna_sequence);
            break;
        }
        

        case 8:
            translateORF(orfs);
            break;
        case 9:
            proteinstatisticsandweight(currentProtein);
            break;
        case 10:
            cout << "Exiting..." << endl;
            return 0;
        default:
            cout << "Invalid Choice!" << endl;
        }
    }
}