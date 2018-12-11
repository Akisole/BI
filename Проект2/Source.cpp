#include <iostream>
#include <string>
using namespace std;

void numsub(string subs, string subf, int &num) {
	int numers = 1;

	for (int i = 0; i <= subf.length() - subs.length(); i++) {

		if (subf[i] == subs[0]) {
			string ssub = subf.substr(i, subs.length());
			if (ssub == subs) {
				numers++;
				i += 4;
			}
		}

	}

	num = numers;
}

void complement(string strin, char &a) {

	if (strin[0] == 'A')	a = 'T';
	else if (strin[0] == 'T')	a = 'A';
	else if (strin[0] == 'G')	a = 'C';
	else if (strin[0] == 'C')	a = 'G';	

}

void codon(string strin, char &a) {

	if (strin == "AAA") a = 'K';
	else if (strin == "AAC") a = 'N';	else if (strin == "AAG") a = 'K';	else if (strin == "AAU") a = 'N';
	else if (strin == "ACA") a = 'T';	else if (strin == "ACC") a = 'T';	else if (strin == "ACG") a = 'T';
	else if (strin == "ACU") a = 'T';	else if (strin == "AGA") a = 'R';	else if (strin == "AGC") a = 'S';
	else if (strin == "AGG") a = 'R';	else if (strin == "AGU") a = 'S';	else if (strin == "AUA") a = 'I';
	else if (strin == "AUC") a = 'I';	else if (strin == "AUG") a = 'M';	else if (strin == "AUU") a = 'I';
	else if (strin == "CAA") a = 'Q';	else if (strin == "CAC") a = 'H';	else if (strin == "CAG") a = 'Q';
	else if (strin == "CAU") a = 'H';	else if (strin == "CCA") a = 'P';	else if (strin == "CCC") a = 'P';
	else if (strin == "CCG") a = 'P';	else if (strin == "CCU") a = 'P';	else if (strin == "CGA") a = 'R';
	else if (strin == "CGC") a = 'R';	else if (strin == "CGG") a = 'R';	else if (strin == "CGU") a = 'R';
	else if (strin == "CUA") a = 'L';	else if (strin == "CUC") a = 'L';	else if (strin == "CUG") a = 'L';
	else if (strin == "CUU") a = 'L';	else if (strin == "GAA") a = 'E';	else if (strin == "GAC") a = 'D';
	else if (strin == "GAG") a = 'E';	else if (strin == "GAU") a = 'D';	else if (strin == "GCA") a = 'A';
	else if (strin == "GCC") a = 'A';	else if (strin == "GCG") a = 'A';	else if (strin == "GCU") a = 'A';
	else if (strin == "GGA") a = 'G';	else if (strin == "GGC") a = 'G';	else if (strin == "GGG") a = 'G';
	else if (strin == "GGU") a = 'G';	else if (strin == "GUA") a = 'V';	else if (strin == "GUC") a = 'V';
	else if (strin == "GUG") a = 'V';	else if (strin == "GUU") a = 'V';	else if (strin == "UAA") a = ' ';
	else if (strin == "UAC") a = 'Y';	else if (strin == "UAG") a = ' ';	else if (strin == "UAU") a = 'Y';
	else if (strin == "UCA") a = 'S';	else if (strin == "UCC") a = 'S';	else if (strin == "UCG") a = 'S';
	else if (strin == "UCU") a = 'S';	else if (strin == "UGA") a = ' ';	else if (strin == "UGC") a = 'C';
	else if (strin == "UGG") a = 'W';	else if (strin == "UGU") a = 'C';	else if (strin == "UUA") a = 'L';
	else if (strin == "UUC") a = 'F';	else if (strin == "UUG") a = 'L';	else if (strin == "UUU") a = 'F';
}

void codon2(string strin1, char &a) {

	string strin;
	for (int i = 0; i < 3; i++) {
		if (strin1[i] == 'T')
			strin += 'U';
		else
			strin += strin1[i];
	}

	if (strin == "AAA") a = 'K';
	else if (strin == "AAC") a = 'N';	else if (strin == "AAG") a = 'K';	else if (strin == "AAU") a = 'N';
	else if (strin == "ACA") a = 'T';	else if (strin == "ACC") a = 'T';	else if (strin == "ACG") a = 'T';
	else if (strin == "ACU") a = 'T';	else if (strin == "AGA") a = 'R';	else if (strin == "AGC") a = 'S';
	else if (strin == "AGG") a = 'R';	else if (strin == "AGU") a = 'S';	else if (strin == "AUA") a = 'I';
	else if (strin == "AUC") a = 'I';	else if (strin == "AUG") a = 'M';	else if (strin == "AUU") a = 'I';
	else if (strin == "CAA") a = 'Q';	else if (strin == "CAC") a = 'H';	else if (strin == "CAG") a = 'Q';
	else if (strin == "CAU") a = 'H';	else if (strin == "CCA") a = 'P';	else if (strin == "CCC") a = 'P';
	else if (strin == "CCG") a = 'P';	else if (strin == "CCU") a = 'P';	else if (strin == "CGA") a = 'R';
	else if (strin == "CGC") a = 'R';	else if (strin == "CGG") a = 'R';	else if (strin == "CGU") a = 'R';
	else if (strin == "CUA") a = 'L';	else if (strin == "CUC") a = 'L';	else if (strin == "CUG") a = 'L';
	else if (strin == "CUU") a = 'L';	else if (strin == "GAA") a = 'E';	else if (strin == "GAC") a = 'D';
	else if (strin == "GAG") a = 'E';	else if (strin == "GAU") a = 'D';	else if (strin == "GCA") a = 'A';
	else if (strin == "GCC") a = 'A';	else if (strin == "GCG") a = 'A';	else if (strin == "GCU") a = 'A';
	else if (strin == "GGA") a = 'G';	else if (strin == "GGC") a = 'G';	else if (strin == "GGG") a = 'G';
	else if (strin == "GGU") a = 'G';	else if (strin == "GUA") a = 'V';	else if (strin == "GUC") a = 'V';
	else if (strin == "GUG") a = 'V';	else if (strin == "GUU") a = 'V';	else if (strin == "UAA") a = ' ';
	else if (strin == "UAC") a = 'Y';	else if (strin == "UAG") a = ' ';	else if (strin == "UAU") a = 'Y';
	else if (strin == "UCA") a = 'S';	else if (strin == "UCC") a = 'S';	else if (strin == "UCG") a = 'S';
	else if (strin == "UCU") a = 'S';	else if (strin == "UGA") a = ' ';	else if (strin == "UGC") a = 'C';
	else if (strin == "UGG") a = 'W';	else if (strin == "UGU") a = 'C';	else if (strin == "UUA") a = 'L';
	else if (strin == "UUC") a = 'F';	else if (strin == "UUG") a = 'L';	else if (strin == "UUU") a = 'F';
}

void reversstr(string strin, string &strout) {
	string revers;
	char a;
	for (int i = 0; i < strin.length(); i++) {
		complement(strin.substr(i, 1), a);
		revers += a;
	}
	for (int i = strin.length() - 1; i >= 0; i--)
		strout += revers[i];
}


int main() {
	// put your code here
	int task;
	
	cout << "Task: ";
	cin >> task;
	
	//--------------------------------//
	//  #1(1.1)	  #2(1.2)	#3(1.3)	  //
	//  #4(2.1)	  #5(2.2)	#6(2.3)	  //
	//  #7(2.4)	  #8(2.5)		      //
	//  #9(3.1)	  #10(3.2)  #11(3.3)  //
	//  #12(4.1)  #13(4.2)			  //
	//--------------------------------//

	while (task != -1) {
		if (task == 1) {
			string str1, str2;
			int numres = 0;

			cin >> str1;
			cin >> str2;

			for (int i = 0; i < str2.length(); i++) {

				if (str2[i] == str1[0]) {
					int ii = i;
					int j = 0;

					while ((ii < str2.length()) && (j < str1.length()) && ((str2[ii] == str1[j]))) {
						ii++;
						j++;
					}
					if (j == str1.length())
						numres++;

				}

			}

			cout << numres;
		}
		if (task == 2) {
			string str;
			int size;

			cin >> str;
			cin >> size;

			string* smas = new string[str.length() - size + 1];


			int n_max = 0;
			int nsub = 0;

			for (int i = 0; i <= str.length() - size; i++) {
				int n;
				string sub = str.substr(i, size);
				if (str.length() - size - i >= 4)
					numsub(sub, str.substr(i + 4, str.length() - size - i), n);
				else if (str.length() - size - i < 4)
					n = 1;

				if (n > n_max) {
					nsub = 0;
					n_max = n;
					smas[nsub] = sub;
					nsub++;
				}
				else if (n == n_max) {
					smas[nsub] = sub;
					nsub++;
				}
			}

			for (int i = 0; i < nsub; i++)
				cout << smas[i] << " ";

			delete[] smas;
		}
		if (task == 3) {
			string str;
			char a;

			cin >> str;

			string* res = new string[str.length()];

			for (int i = 0; i < str.length(); i++) {
				complement(str.substr(i, 1), a);
				res[i] = a;
			}

			for (int i = str.length() - 1; i >= 0; i--)
				cout << res[i];

			delete[] res;
		}
		if (task == 4) {

			string str, res;

			cin >> str;

			for (int i = 0; i <= str.length() - 3; i = i + 3) {
				char a;
				codon(str.substr(i, 3), a);
				res += a;
			}

			cout << res;
		}
		if (task == 5) {

			string str, strf;

			cin >> str;
			cin >> strf;

			string* res = new string[str.length()];
			int num_r = 0;

			for (int i = 0; i <= str.length() - 6; i++) {
				char ach, bch;
				string ast = str.substr(i, 3);
				string bst = str.substr(i + 3, 3);
				string strtmp = ast + bst;

				codon2(ast, ach);
				codon2(bst, bch);

				if (strf[0] == ach && strf[1] == bch) {
					res[num_r] = ast + bst;
					num_r++;
				}

				string revers;
				reversstr(strtmp, revers);
				ast = revers.substr(0, 3);
				bst = revers.substr(3, 3);
				codon2(ast, ach);
				codon2(bst, bch);
				if (strf[0] == ach && strf[1] == bch) {
					res[num_r] = strtmp;
					num_r++;
				}
			}


			for (int i = 0; i < num_r; i++)
				cout << res[i] << endl;


		}
		if (task == 6) {
			int num;
			cin >> num;
			cout << (num - 1)*num;
		}

		cout << endl << "Task: ";
		cin >> task;
	}
	


	return 0;
}