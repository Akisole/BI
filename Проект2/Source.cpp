#include <iostream>
#include <string>
#include <vector>
#include <map> 
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

void alfamin(string strin1, int &a) {	
	if (strin1 == "G") a = 57;			else if (strin1 == "A") a = 71;		else if (strin1 == "S") a = 87; 
	else if (strin1 == "P") a = 97;		else if (strin1 == "V") a = 99;		else if (strin1 == "T") a = 101;
	else if (strin1 == "C") a = 103;	else if (strin1 == "I") a = 113;	else if (strin1 == "L") a = 113;
	else if (strin1 == "N") a = 114;	else if (strin1 == "D") a = 115;	else if (strin1 == "K") a = 128;
	else if (strin1 == "Q") a = 128;	else if (strin1 == "E") a = 129;	else if (strin1 == "M") a = 131;
	else if (strin1 == "H") a = 137;	else if (strin1 == "F") a = 147;	else if (strin1 == "R") a = 156;
	else if (strin1 == "Y") a = 163;	else if (strin1 == "W") a = 186;
}
string alfaminum(int i) {
	if (i == 1) return "G";
	else if (i == 2) return "A";
	else if (i == 3) return "S";
	else if (i == 4) return "P";
	else if (i == 5) return "V";
	else if (i == 6) return "T";
	else if (i == 7) return "C";
	else if (i == 8) return "I";
	/*else if (strin1 == "L") a = 113;*/
	else if (i == 9) return "N";
	else if (i == 10) return "D";
	else if (i == 11) return "K";
	/*else if (strin1 == "Q") a = 128;*/
	else if (i == 12) return "E";
	else if (i == 13) return "M";
	else if (i == 14) return "H";
	else if (i == 15) return "F";
	else if (i == 16) return "R";
	else if (i == 17) return "Y";
	else if (i == 18) return "W";
}

string alfamin(int a) {
	if (a == 57) return "G";
	else if (a == 71) return "A";
	else if (a == 87) return "S";
	else if (a == 97) return "P";
	else if (a == 99) return "V";		
	else if (a == 101) return "T";
	else if (a == 103) return "C";	
	else if (a == 113) return "I";	
	/*else if (strin1 == "L") a = 113;*/
	else if (a == 114) return "N";	
	else if (a == 115) return "D";	
	else if (a == 128) return "K";
	/*else if (strin1 == "Q") a = 128;*/	
	else if (a == 129) return "E";
	else if (a == 131) return "M";
	else if (a == 137) return "H";
	else if (a == 147) return "F";
	else if (a == 156) return "R";
	else if (a == 163) return "Y";
	else if (a == 186) return "W";
}

void sort(int *mas, int len) {
	int i, j, b;
	for (i = 0; i<len; i++)
		for (j = 0; j<len - 1 - i; j++)
		{
			if (mas[j]>mas[j + 1])
			{
				b = mas[j];
				mas[j] = mas[j + 1];
				mas[j + 1] = b;
			}
		}
}

void sort(string *mas, int len) {
	int i, j;
	string b;
	for (i = 0; i<len; i++)
		for (j = 0; j<len - 1 - i; j++)
		{
			if (mas[j]>mas[j + 1])
			{
				b = mas[j];
				mas[j] = mas[j + 1];
				mas[j + 1] = b;
			}
		}
}

int* delrepitesort(int *mas, int len, int &numel) {
	int* res;
	int num=1, a=mas[0], k;

	for (int i=1; i<len;i++)
		if (mas[i] > a) {
			num++;
			a = mas[i];
		}
	res = new int[num];

	a = mas[0];
	k = 0;
	res[k++] = a;
	for (int i = 1; i<len; i++)
		if (mas[i] > a) {
			a = mas[i];
			res[k++] = a;
		}
	numel = num;
	return res;
}

int mass_peptide(string strin) {
	int mas = 0, a;
	for (int i = 0; i < strin.length(); i++) {
		alfamin(strin.substr(i, 1), a);
		mas += a;
	}
	return mas;
}

int* StrToInt(string &str, int &size) {

	string strk;
	int i, p;
	int *mass;
	getline(cin, str);

	size = 0;
	for (i = 0; i < str.length(); i = p + 1) {
		strk = str.substr(i, str.length() - i);
		p = strk.find(' ');
		if (p > 0 && p < str.length()) {
			p = p + i;
			size++;
		}
		else
			break;
	}
	size++;

	mass = new int[size];

	int k = 0;
	for (i = 0; i < str.length(); i = p + i + 1) {
		strk = str.substr(i, str.length() - i);
		p = strk.find(' ');
		if (p != -1) {
			mass[k] = stoi(str.substr(i, p));
			k++;
		}
		else
			break;
	}
	mass[k] = stoi(strk);

	return mass;
}

int* Cyclospectrum(string strin, int &size) {
	string  strk;
	int *masami, *circl;
	int circllen, i, j, k, a;

	size = strin.length()*(strin.length() - 1) + 2;

	masami = new int[size];


	circllen = strin.length() + strin.length() - 2;
	if (circllen>0)
		circl = new int[circllen];


	for (i = 0; i < strin.length(); i++) {
		alfamin(strin.substr(i, 1), a);
		masami[i] = a;
		if (circllen>0)
			circl[i] = a;
	}

	j = 0;
	if (circllen>0)
		for (i; i < circllen; i++)
			circl[i] = masami[j++];


	k = strin.length();
	for (j = 2; j < strin.length(); j++) {
		for (i = 0; i < strin.length(); i++) {

			a = 0;
			if (circllen>0)
				for (int ii = i; ii < (i + j); ii++)
					a += circl[ii];
			masami[k] = a;
			k++;
		}
	}

	masami[k] = 0;
	a = 0;
	for (i = 0; i < strin.length(); i++)
		a += masami[i];
	masami[k + 1] = a;

	sort(masami, size);

	if (circllen>0)
		delete[]circl;
	return masami;
}

int* TheorLinerspectrum(string strin, int &size) {

	string  strk;
	int *masami, *circl;
	int circllen, i, j, k, a;

	size = 1;
	for (i = 1; i <= strin.length(); i++) 
		size += i;
	
	masami = new int[size];

	for (i = 0; i < strin.length(); i++) {
		alfamin(strin.substr(i, 1), a);
		masami[i] = a;
	}

	k = i;
	for (j = 2; j <= strin.length(); j++) {
		
		for (int i = 0; i <= strin.length() - j; i++) {
			a = 0;
			for (int ii = i; ii < (i + j); ii++)
				a += masami[ii];
			masami[k] = a;
			k++;
		}
	}

	masami[k] = 0;
	sort(masami, size);

	return masami;
}

bool equalCyclSpPept_Spectrum(int* SpectPept, int s_SP, int* Spect, int s_S) {
	if (s_S != s_SP)
		return false;
	for (int i = 0; i < s_S; i++)
		if (Spect[i] != SpectPept[i])
			return false;
	return true;
}

bool consistentPept_Spect(int* SpectPept, int s_SP, int* Spect, int s_S) {
	int score, i, j;

	score = i = j = 0;
	while (i < s_SP && j < s_S) {
		if (SpectPept[i] == Spect[j]) {
			i++;
			j++;
			score++;
		}
		else if (SpectPept[i] > Spect[j])
			j++;
		else
			i++;
	}

	if (score == s_SP)
		return true;

	return false;
}

int removePept(string* peptides, int numofpept) {
	int size = numofpept;
	for (int i = 0; i < size; i++)
		if (peptides[i] == "-1") {
			peptides[i] = peptides[size - 1];
			size--;
			i--;
		}

	return size;
}

void expend(string* peptides, int &numofpept) {
	int i, j, size;
	size = numofpept;

	for (i = 0; i < size; i++) 
		for (j = 1; j <= 18; j++) 
			peptides[numofpept++] = peptides[i] + alfaminum(j);		

	for (i = 0; i < numofpept - size; i++)
		peptides[i] = peptides[i + size];

	numofpept -= size;
}

void answerSequencing(string peptides) {
	string ans="";
	int p;
	for (int i = 0; i < peptides.length(); i++) {
		if (peptides.substr(i, 1) == " ") {
			ans.pop_back();
			ans += " ";
		}
		else {
			int a;
			alfamin(peptides.substr(i, 1), a);
			ans += to_string(a) + '-';
		}
	}

	if (ans[ans.length()-1] == '-')
		ans.pop_back();
	cout << ans;
}

int Score(int* PeptidCyclo, int s_PC, int* Spectrum, int s_S) {
	int i=0, j=0, score = 0;

	while (i < s_PC && j < s_S) {
		if (PeptidCyclo[i] == Spectrum[j]) {
			i++;
			j++;
			score++;
		}
		else if (PeptidCyclo[i] > Spectrum [j])
			j++;
		else
			i++;
	}
	return score;
}

int LinearScore(int* PeptidLiner, int s_PC, int* Spectrum, int s_S) {
	int i=0, j=0, score = 0;

	while (i < s_PC && j < s_S) {
		if (PeptidLiner[i] == Spectrum[j]) {
			i++;
			j++;
			score++;
		}
		else if (PeptidLiner[i] > Spectrum[j])
			j++;
		else
			i++;
	}
	return score;
}

void Tream(string* leadbord, int &numofpept, int* Spectrum, int s_S, int N) {
	if (numofpept == 0) return;

	int max_scor=0, prev_scor = 0, num_maxscor = 0, a, i, j, k, in, n = N, size_wr=0;
	int* scorpepts = new int[numofpept];
	int* scorpepts_sort = new int[numofpept];
	int* scorsort_wrepit;
	

	string* tmp;

	for (i = 0; i < numofpept; i++) {
		scorpepts_sort[i] = LinearScore(TheorLinerspectrum(leadbord[i], a), a, Spectrum, s_S);
		scorpepts[i] = LinearScore(TheorLinerspectrum(leadbord[i], a), a, Spectrum, s_S);
		if (scorpepts[i] > max_scor) {
			prev_scor = max_scor;
			max_scor = scorpepts[i];
		}
	}

	sort(scorpepts_sort, numofpept);
	scorsort_wrepit = delrepitesort(scorpepts_sort, numofpept, size_wr);

	for (i = 0; i < numofpept; i++) {
		if (scorpepts[i] == max_scor)
			num_maxscor++;
	}

	if (num_maxscor >= N) {
		tmp = new string[num_maxscor];
		for (i = 0, in = 0; i < numofpept && in < num_maxscor; i++)
			if (scorpepts[i] == max_scor) {
				tmp[in++] = leadbord[i];
			}
		for (i = 0; i < num_maxscor; i++)
			leadbord[i] = tmp[i];
		numofpept = num_maxscor;
	}
	else {
		tmp = new string[N];
		for (i = 0, in = 0; i < numofpept && in < num_maxscor; i++)
			if (scorpepts[i] == max_scor) {
				tmp[in++] = leadbord[i];
			}

		for (j = size_wr - 2; j >= 0; j++)
			for (i = 0; i < numofpept, in<N; i++)
				if (scorpepts[i] == scorsort_wrepit[j])
					tmp[in++] = leadbord[i];

		for (i = 0; i < N; i++)
			leadbord[i] = tmp[i];
		numofpept = N;
	}
	
}
//12
bool mistakes_d(string pat1, string pat2, int k, int d) {
	int score = 0;
	for (int i = 0; i < k; i++) {
		if (pat1[i] != pat2[i])
			score++;
	}
	if (score == d || score < d)
		return true;
	return false;

}

bool checkstrdna_d(string str, string pat, int k, int d) {
	
	for (int i = 0; i < str.length() - k + 1; i++) {
		string tmp = str.substr(i, k);
		if (mistakes_d(pat, tmp, k, d))
			return true;
	}
	return false;
}

string putmistakes(string pat, int d, int i, int l) {
	int ll;
	string res = pat;
	char a = pat[i];
	if (a == 'A') ll = l + 2;
	else if (a == 'C') ll = l + 3;
	else if (a == 'T') ll = l + 1;
	else if (a == 'G') ll = l + 4;

	if (d == 1) {
		if (ll == 1) {
			res[i] = 'A';
			return res;
		}
		if (ll == 2) {
			res[i] = 'C';
			return res;
		}
		if (ll == 3) {
			res[i] = 'G';
			return res;
		}
		if (ll == 4) {
			res[i] = 'T';
			return res;
		}
		if (ll == 5) {
			res[i] = 'A';
			return res;
		}
		if (ll == 6) {
			res[i] = 'C';
			return res;
		}
	}

}
//13
vector<string> abc{ "A", "C", "G", "T" };
void gen(const vector<string> &alf, string s, size_t n, string* mas, int &mcount) {
	for (size_t i = 0; i < alf.size(); ++i) {
		if (n > 1) 
			gen(alf, s + alf[i], n - 1,mas,mcount);
		else 
			mas[mcount++] = s + alf[i];
	}
}

int HammingDist(string pat, string subpat)
{
	int count = 0;
	int length = pat.length();

	for (int i = 0; i < length; i++){
		if (pat[i] != subpat[i])
			count++;
	}
	return count;
}

int Dpat_text(string pat, string dnai, int k) {
	int l = dnai.length();
	int d = k, tmpd;

	for (int i = 0; i < l - k + 1; i++) {
		tmpd = HammingDist(pat, dnai.substr(i, k));
		if (tmpd < d)
			d = tmpd;
	}
	return d;
}

int Dpat_dna(string pat, string* dna, int k) {
	int sum = 0;
	for (int i = 0; i < 10; i++) 
		sum += Dpat_text(pat, dna[i], k);
	return sum;
}
string* delrepit(string str, int k, int &size) {
	string strk, strk2, strpat, strres="", *resmas;
	int i, j, flag, massize = 0;

	resmas = new string[str.length() / k];

	strk = str;
	for (int i = 0; i < str.length() - k + 1; i=i+k+1) {
		strk = str.substr(i, k);
		flag = 0;
		for (j = 0; j < massize; j++) {
			if (resmas[j] == strk)
				flag++;
		}
		if (!flag) {
			resmas[massize++] = strk;
		}
	}

	size = massize;
	return resmas;
}



int main() {
	// put your code here
	int task;

	cout << "Task: ";
	cin >> task;

	//--------------------------------//
	//  #1(1.1)	  #2(1.2)	#3(1.3)	  //
	//  #4(2.1)	  #5(2.2)	#6(2.3)	  //
	//--------------------------------//
	//  #7(2.4)	  #8(2.5)		      //
	//  #9(3.1)  #10(3.2)  #11(3.3)  //
	//  #12(4.1) #13(4.2)		 //
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
		if (task == 7) {
			string str1, strk;
			int *masami, *circl;
			int size, circllen, i, j, a;

			cin >> str1;
			size = str1.length()*(str1.length() - 1) + 2;

			masami = new int[size];

			circllen = str1.length() + str1.length() - 2;
			circl = new int[circllen];


			for (i = 0; i < str1.length(); i++) {
				alfamin(str1.substr(i, 1), a);
				masami[i] = a;
				circl[i] = a;
			}

			j = 0;
			for (i; i < circllen; i++)
				circl[i] = masami[j++];


			int k = str1.length();
			for (j = 2; j < str1.length(); j++) {
				for (i = 0; i < str1.length(); i++) {

					a = 0;
					for (int ii = i; ii < (i + j); ii++)
						a += circl[ii];
					masami[k] = a;
					k++;
				}
			}

			masami[k] = 0;
			a = 0;
			for (i = 0; i < str1.length(); i++)
				a += masami[i];
			masami[k + 1] = a;

			sort(masami, size);
			for (i = 0; i < size; i++)
				cout << masami[i] << " ";
		}
		if (task == 8) {

			int aminoacid_masses[18] = { 57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186 };

			map<int, int64_t> res = { { 0,1 } };
			int input;
			cin >> input;

			for (int i = 57; i <= input; i++)
			{
				res[i] = 0;

				for (int j = 0; j < 18; j++)
				{
					if (res.find(i - aminoacid_masses[j]) != res.end())
					{
						res[i] += res[i - aminoacid_masses[j]];
					}
				}
			}
			cout << res[input];
		}
		if (task == 9) {
			string str1, *strvarians, strresultpepts = "";
			int *spectr;
			int size = 0, parentmass_spctr, size_outpept, numvarianspept, currvarianse;
			int i;

			cin.get();
			spectr = StrToInt(str1, size);
			parentmass_spctr = spectr[size - 1];

			for (i = 0; i < size; i++) {
				if (i * (i + 1) == size - 2) {
					size_outpept = i + 1;
					break;
				}
			}

			numvarianspept = pow(18, size_outpept);

			strvarians = new string[numvarianspept];

			currvarianse = 0;

			for (i = 1; i <= 18; i++) {
				strvarians[i - 1] = alfaminum(i);
				currvarianse++;
			}

			int flag = 0;
			while (currvarianse != 0) {

				if (flag)
					expend(strvarians, currvarianse);

				for (i = 0; i < currvarianse; i++) {
					int a;
					if (mass_peptide(strvarians[i]) == parentmass_spctr) {

						if (equalCyclSpPept_Spectrum(Cyclospectrum(strvarians[i], a), a, spectr, size))
							strresultpepts += strvarians[i] + ' ';

						strvarians[i] = "-1";
					}
					else if (!consistentPept_Spect(Cyclospectrum(strvarians[i], a), a, spectr, size))
						strvarians[i] = "-1";
				}

				currvarianse = removePept(strvarians, currvarianse);

				if (!flag)	flag = 1;
			}

			answerSequencing(strresultpepts);
			
		}
		if (task == 10) {
			string str1, str2;
			int *masspectr, *peptidspectr;
			int sizesp, sizepeptsp, i, j, score;

			cin >> str1;
			cin.get();
			masspectr = StrToInt(str2, sizesp);
			peptidspectr = Cyclospectrum(str1, sizepeptsp);

			score = i = j = 0;
			while (i < sizesp && j < sizepeptsp) {
				if (masspectr[i] == peptidspectr[j]) {
					i++;
					j++;
					score++;
				}
				else if (masspectr[i] > peptidspectr[j])
					j++;
				else
					i++;
			}

			cout << score;

			delete[] masspectr;
			delete[] peptidspectr;

		}
		if (task == 11) {
			string str1, *strvarians, resultpepts = "";
			int *spectr;
			int size = 0, parentmass_spctr, size_outpept, numvarianspept, currvarianse;
			int i, N;

			cin >> N;

			cin.get();
			spectr = StrToInt(str1, size);
			parentmass_spctr = spectr[size - 1];

			size_outpept = N;

			numvarianspept = pow(18, 4);
			strvarians = new string[numvarianspept];
			currvarianse = 0;

			for (i = 1; i <= 18; i++) {
				strvarians[i - 1] = alfaminum(i);
				currvarianse++;
			}

			int flag = 0;
			while (currvarianse != 0) {
				expend(strvarians, currvarianse);


				for (i = 0; i < currvarianse; i++) {
					int a, b;
					if (mass_peptide(strvarians[i]) == parentmass_spctr) {
						if (LinearScore(TheorLinerspectrum(strvarians[i], a), a, spectr, size) >
							LinearScore(TheorLinerspectrum(resultpepts, b), b, spectr, size)) {

							resultpepts = strvarians[i];
						}
					}
					else if (mass_peptide(strvarians[i]) > parentmass_spctr)
						strvarians[i] = "-1";	
				}


				currvarianse = removePept(strvarians, currvarianse);

				Tream(strvarians, currvarianse, spectr, size, N);

				if (!flag)	flag = 1;
			}

			answerSequencing(resultpepts);

		}
		if (task == 12) {
			string *dna, patern, tmp, tmppat, result="", *resmas;
			int k, d, i, jd, ii,j;
			int *mastmp, size, checknum, mistnam, rmsize;
			bool check;
			char a;

			cin.get();
			mastmp = StrToInt(tmp, size);
			dna = new string[4];
			for (i = 0; i < 4; i++)
				cin >> dna[i];

			k = mastmp[0];
			d = mastmp[1];

			for (i = 0; i < 4; i++) {
				for (jd = 0; jd < dna[i].length() - k + 1; jd++) {
					tmppat = dna[i].substr(jd, k);
					
					
					for (j = 0; j < k; j++) {
						for (mistnam = 0; mistnam < 3; mistnam++) {

							string patm;
							patm = putmistakes(tmppat, d, j, mistnam);
							checknum = 0;
							for (ii = 0; ii < 4; ii++) {
								if (checkstrdna_d(dna[ii], patm, k, d))
									checknum++;
							}
							if (checknum == 4)
								result += patm + ' ';
						}
					}
				}
			}

			result.pop_back();

			resmas = delrepit(result, k, rmsize);
			
			sort(resmas, rmsize);

			for (i = 0; i < rmsize; i++)
				cout << resmas[i] << " ";

		}
		if (task == 13) {
			int i, distanse = 1000, tmpd;
			int k, mascount = 0;
			string *mass, s, *dna, median;

			cin >> k;
			dna = new string[10];
			for (i = 0; i < 10; i++)
				cin >> dna[i];

			mass = new string[pow(4, k)];
			gen(abc, s, k,mass,mascount);

			for (i = 0; i < mascount; i++) {
				tmpd = Dpat_dna(mass[i], dna, k);
				if (distanse > tmpd) {
					distanse = tmpd;
					median = mass[i];
				}
			}
			cout << median;

		}


		cout << endl << "Task: ";
		cin >> task;
	}



	return 0;
}
