#include <iostream>
#include <string>
using namespace std;

void complement(string strin, char &a) {

	if (strin[0] == 'A')	a = 'T';
	else if (strin[0] == 'T')	a = 'A';
	else if (strin[0] == 'G')	a = 'C';
	else if (strin[0] == 'C')	a = 'G';

}
void numsub(string subs, string subf, int &num) {
	int numers = 1;

	if (subf.length()<subs.length())
		return;

	for (int i = 0; i <= subf.length() - subs.length(); i++) {

		if (subf[i] == subs[0]) {
			string ssub = subf.substr(i, subs.length());
			if (ssub == subs) {
				numers++;
				//				i += subs.length();
			}
		}
	}
	num = numers;
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
string codon2str(string str) {
	string res = "";
	char a;
	for (int i = 0; i < str.length(); i += 3) {
		codon2(str.substr(i, 3), a);
		res += a;
	}
	return res;
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
void alfamin(char strin1, int &a) {
	if (strin1 == 'G') a = 57;			else if (strin1 == 'A') a = 71;		else if (strin1 == 'S') a = 87;
	else if (strin1 == 'P') a = 97;		else if (strin1 == 'V') a = 99;		else if (strin1 == 'T') a = 101;
	else if (strin1 == 'C') a = 103;	else if (strin1 == 'I') a = 113;	else if (strin1 == 'N') a = 114;
	else if (strin1 == 'D') a = 115;	else if (strin1 == 'K') a = 128;	else if (strin1 == 'E') a = 129;
	else if (strin1 == 'M') a = 131;	else if (strin1 == 'H') a = 137;	else if (strin1 == 'F') a = 147;
	else if (strin1 == 'R') a = 156;	else if (strin1 == 'Y') a = 163;	else if (strin1 == 'W') a = 186;
}
string alfamin(int a) {
	if (a == 57) return "G"; 	else if (a == 71) return "A"; 	else if (a == 87) return "S";
	else if (a == 97) return "P";    else if (a == 99) return "V";	else if (a == 101) return "T";
	else if (a == 103) return "C";	else if (a == 113) return "I";	else if (a == 114) return "N";
	else if (a == 115) return "D";	else if (a == 128) return "K";   else if (a == 129) return "E";
	else if (a == 131) return "M";	else if (a == 137) return "H";	else if (a == 147) return "F";
	else if (a == 156) return "R";	else if (a == 163) return "Y";	else if (a == 186) return "W";
}
string alfaminum(int i) {
	if (i == 1) return "G"; 	else if (i == 2) return "A"; else if (i == 3) return "S";
	else if (i == 4) return "P";	else if (i == 5) return "V";	else if (i == 6) return "T";
	else if (i == 7) return "C";	else if (i == 8) return "I";	else if (i == 9) return "N";
	else if (i == 10) return "D";	else if (i == 11) return "K";	else if (i == 12) return "E";
	else if (i == 13) return "M";	else if (i == 14) return "H";	else if (i == 15) return "F";
	else if (i == 16) return "R";	else if (i == 17) return "Y";	else if (i == 18) return "W";
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
string delrepitesort(string *mas, int len) {
	string res = "";
	int num = 1, k = 0, i;
	string a = mas[0];

	for (i = 0; i < len - 1; i++) {
		a = mas[i];
		int flag = 0;
		for (int j = i + 1; j < len; j++) {
			if (a == mas[j]) {
				flag = 1;
				break;
			}
		}
		if (!flag) {
			res += mas[i] + ' ';
		}
	}
	res += mas[i];
	return res;
}

int* StrToInt(string &str, int &size) {

	string strk;
	int i, p;
	int *mass;
	getline(cin, str);
/*
	if (str[str.length() - 1] == ' ') {
		while (str[str.length() - 1] == ' ') {
			str.pop_back();
		}
	}
	if (str[0] == ' ') {
		while (str[0] == ' ') {
			str = str.substr(1, str.length() - 1);
		}
	}
*/
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
void answerSequencing(string *peptides, int size) {
	string ans = "";
	int mas, size_wr;
	string *scorsort_wrepit;

	for (int i = 0; i < size; i++) {
		mas = 0;
		ans = "";
		for (int j = 0; j < peptides[i].length(); j++) {
			alfamin(peptides[i][j], mas);
			ans += to_string(mas) + '-';
		}
		ans.pop_back();
		peptides[i] = ans;
	}
	sort(peptides, size);

	cout << delrepitesort(peptides, size);

}
void answerSequencing(string peptides, int varians) {
	/*string ans="";
	string *ansmas;
	ansmas = new string[varians];
	int p, massize=0;
	for (int i = 0; i < peptides.length(); i++) {
	if (peptides.substr(i, 1) == " ") {
	ans.pop_back();
	ansmas[massize++] = ans;
	//			ans += " ";
	ans = "";
	}
	else {
	int a;
	alfamin(peptides.substr(i, 1), a);
	ans += to_string(a) + '-';

	}
	}
	sort(ansmas, massize);
	for (int i=0; i<varians; i++)
	cout << ansmas[i] << " ";
	delete[] ansmas;*/

	//	if (ans[ans.length()-1] == '-')
	//		ans.pop_back();
	//	cout << ans;
}
int mass_peptide(string strin) {
	int mas = 0, a;
	for (int i = 0; i < strin.length(); i++) {
		alfamin(strin.substr(i, 1), a);
		mas += a;
	}
	return mas;
}

int* TheorLinerspectrum(string strin, int &size) {

	string  strk;
	int *masami;
	int i, j, k, a;

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

int* delrepitesort(int *mas, int len, int &numel) {
	int* res;
	int num = 1, a = mas[0], k;

	for (int i = 1; i<len; i++)
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
int LinearScore(int* PeptidLiner, int s_PC, int* Spectrum, int s_S) {
	int i = 0, j = 0, score = 0;
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

	int max_scor = 0, prev_scor = 0, num_maxscor = 0, a, i, j, k, in, n = N, size_wr = 0;
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
void answerSequencing(string peptides) {
	string ans = "";
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
	if (ans[ans.length() - 1] == '-')
		ans.pop_back();
	cout << ans;
}

int main() {
	// put your code here
	int task;

	cout << "Task (2/5/6/9/10/11): ";
	cin >> task;

	while (task != -1) {
		if (task == 2) {
			cout << "Task 1.2" << endl;
			string str, strsub;
			long int k;

			cin >> str;
			cin >> k;

			string* smas = new string[str.length()];

			int n_max = 0;
			int nsub = 0;

			for (int i = 0; i <= str.length() - k; i++) {
				int n = 0;
				string sub = str.substr(i, k);
				strsub = str.substr(i, str.length() - 1 - i);

				numsub(sub, strsub, n);

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
			sort(smas, nsub);

			for (int i = 0; i < nsub; i++)
				cout << smas[i] << " ";

			delete[] smas;

		}
		if (task == 5) {
			cout << "Task 2.2" << endl;
			string str, amin;

			cin >> str;
			cin >> amin;

			string* res = new string[str.length()];
			int num_r = 0, size_rpep = amin.length() * 3;

			for (int i = 0; i <= str.length() - size_rpep; i++) {

				string strtmp = str.substr(i, size_rpep);
				string tmpamin = codon2str(strtmp);
				if (tmpamin == amin) {
					res[num_r] = strtmp;
					num_r++;
				}

				string revers;
				reversstr(strtmp, revers);
				tmpamin = codon2str(revers);
				if (tmpamin == amin) {
					res[num_r] = strtmp;
					num_r++;
				}

			}

			sort(res, num_r);
			for (int i = 0; i < num_r; i++)
				cout << res[i] << endl;
		}
		if (task == 6) {
			cout << "Task 2.3" << endl;
			long int num;
			cin >> num;
			num = num*(num - 1);
			cout << num;
		}
		if (task == 9) {
			cout << "Task 3.1" << endl;
			string str1, *strvarians, *strresultpeptsmas, strresultpepts = "";
			int *spectr;
			int size = 0, parentmass_spctr, size_outpept, numvarianspept, currvarianse;
			int i, varians = 0;

			cin.get();
			spectr = StrToInt(str1, size);
			sort(spectr, size);
			parentmass_spctr = spectr[size - 1];

			numvarianspept = pow(18, 3);

			strvarians = new string[numvarianspept];
			strresultpeptsmas = new string[numvarianspept];
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
					int a, b;
					if (mass_peptide(strvarians[i]) == parentmass_spctr) {

						int *tmpcyclo = Cyclospectrum(strvarians[i], a);

						if (equalCyclSpPept_Spectrum(tmpcyclo, a, spectr, size)) {
							strresultpepts += strvarians[i] + ' ';
							strresultpeptsmas[varians] = strvarians[i];
							varians++;
						}

						delete[] tmpcyclo;

						strvarians[i] = "-1";
					}
					else {
						int *tmplinar = TheorLinerspectrum(strvarians[i], b);
						if (!consistentPept_Spect(tmplinar, b, spectr, size)) {
							strvarians[i] = "-1";
						}
						delete[] tmplinar;
					}
				}

				currvarianse = removePept(strvarians, currvarianse);

				if (!flag)	flag = 1;
			}

			answerSequencing(strresultpepts, varians);
			//	   cout << endl;
			answerSequencing(strresultpeptsmas, varians);

			delete[] spectr;
			delete[] strvarians;
		}
		if (task == 10) {
			cout << "Task 3.2" << endl;
			string str1, str2;
			int *masspectr, *peptidspectr;
			int sizesp, sizepeptsp, i, j, score;

			cin >> str1;
			cin.get();
			masspectr = StrToInt(str2, sizesp);
			peptidspectr = Cyclospectrum(str1, sizepeptsp);


			for (int i = 0; i < sizepeptsp; i++)
				cout << peptidspectr[i] << " ";
			cout << endl;

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
			cout << "Task 3.3" << endl;
			string str1, *strvarians, resultpepts = "";
			int *spectr;
			int size = 0, parentmass_spctr, size_outpept, currvarianse;
			long int numvarianspept;
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
		cout << endl << "Task: ";
		cin >> task;
	}
	return 0;


}