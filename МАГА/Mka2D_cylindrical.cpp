#include "Mka2D_cylindrical.h"

Mka2D_cylindrical::Mka2D_cylindrical() :Mka2D_cylindrical("resources2D/", false, true, false, false)
{
}

Mka2D_cylindrical::Mka2D_cylindrical(string filePrefix, bool ku1_left, bool ku1_right, bool ku1_up, bool ku1_down) :
	filePrefix(filePrefix), ku1_left(ku1_left), ku1_right(ku1_right), ku1_up(ku1_up), ku1_down(ku1_down)
{
	testFile.open(filePrefix + "log.txt");
}

void Mka2D_cylindrical::startFull2dProcess() {
	GenerateNet();
	directSolveStraightTask();
}

Mka2D_cylindrical::~Mka2D_cylindrical()
{
}


void GenerateNetLikeTelma(set<double> &mas, ifstream &fileNet)
{
	double firstElement, startPosition, endPosition, position, lastPosition, stepReal;;
	int i, numberOfInterval;
	fileNet >> firstElement >> numberOfInterval;
	double*intervals = new double[numberOfInterval + 1];
	double*sizeOfSteps = new double[numberOfInterval];
	double*mnojiteli = new double[numberOfInterval];
	int*napravlenie = new int[numberOfInterval];
	intervals[0] = firstElement;
	for (i = 1; i < numberOfInterval + 1; i++)
	{
		fileNet >> intervals[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> sizeOfSteps[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> mnojiteli[i];
		if (mnojiteli[i] < 1)
		{
			cout << "Ошибка описания сетки по оси. Коэфициент разрядки должен быть больше или равен 1";
			system("pause");
			exit(1);
		}
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> napravlenie[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		if (napravlenie[i] == -1)
		{
			startPosition = intervals[i + 1];
			endPosition = intervals[i];
		}
		else
		{
			startPosition = intervals[i];
			endPosition = intervals[i + 1];
		}
		stepReal = 0;
		lastPosition = startPosition;
		position = startPosition + napravlenie[i] * sizeOfSteps[i];
		while (position*napravlenie[i] < endPosition*napravlenie[i])
		{
			mas.insert(position);
			stepReal = fabs(lastPosition - position)*mnojiteli[i];
			lastPosition = position;
			position += napravlenie[i] * stepReal;
		}
		//if (fabs(lastPosition - endPosition) < sizeOfSteps[i] && lastPosition != startPosition)
		//{
		//	mas.erase(mas.find(lastPosition));
		//}
	}
}

int Mka2D_cylindrical::FindAreaNumber(int nodes[])
{
	int i;
	for (i = 0; i < sreda.size(); i++)
	{
		if (rz[nodes[0]].r >= sreda[i].r1 && rz[nodes[1]].r <= sreda[i].r2 && rz[nodes[0]].z >= sreda[i].z1 && rz[nodes[2]].z <= sreda[i].z2)
			return i;
	}
	cout << "Ошибка в FindAreaNumber: не найдена подобласть." << endl;
	throw new exception("Ошибка в FindAreaNumber: не найдена подобласть.");
	return -1;
}

void Mka2D_cylindrical::GenerateNet()
{
	int i, j, k, numberFields, m;
	int locateSourceX = 0, locateSourceY = 0;
	field tempField;
	set<double> xTemp, yTemp;
	set<double>::const_iterator it;
	ifstream inpSreda(filePrefix + "sreda2D.txt");
	inpSreda >> numberFields;
	double areaOfSreda = 0;
	for (i = 0; i < numberFields; i++)
	{
		inpSreda >> tempField.r1 >> tempField.r2 >> tempField.z1 >> tempField.z2 >> tempField.sigma >> tempField.gamma;
		sreda.push_back(tempField);
		xTemp.insert(tempField.r1);
		xTemp.insert(tempField.r2);
		yTemp.insert(tempField.z1);
		yTemp.insert(tempField.z2);
		areaOfSreda += (tempField.r2 - tempField.r1)*(tempField.z2 - tempField.z1);
	}
	GenerateNetLikeTelma(xTemp, inpSreda);
	GenerateNetLikeTelma(yTemp, inpSreda);
	set<double>::reverse_iterator rit;
	it = xTemp.begin();
	r0 = *it;
	rit = xTemp.rbegin();
	r1 = *rit;

	it = yTemp.begin();
	z0 = *it;
	rit = yTemp.rbegin();
	z1 = *rit;

	for (set<double>::const_iterator it = xTemp.begin(); it != xTemp.end(); it++, i++)
	{
		R.push_back(*it);
	}
	i = 0;
	for (set<double>::const_iterator it = yTemp.begin(); it != yTemp.end(); it++, i++)
	{
		Z.push_back(*it);
	}

	int iKE = 0, iPoint = 0;
	for (size_t iY = 0; iY < Z.size(); iY++)
	{
		for (size_t iX = 0; iX < R.size(); iX++, iPoint++)
		{
			rz.push_back(Point_cylindrical(R[iX], Z[iY]));
		}
	}
	nPoints = rz.size();
	for (size_t iY = 0; iY < Z.size() - 1; iY++)
	{
		for (size_t iX = 0; iX < R.size() - 1; iX++)
		{
			nvtr tmp;
			tmp.uzel[0] = indexRZ(R[iX], Z[iY]);

			tmp.uzel[1] = tmp.uzel[0] + 1;

			tmp.uzel[2] = indexRZ(R[iX], Z[iY + 1]);

			tmp.uzel[3] = tmp.uzel[2] + 1;
			tmp.numberField = FindAreaNumber(tmp.uzel);
			KE.push_back(tmp);
		}
	}

	ofstream fileXY(filePrefix + "rz.txt");
	ofstream fileNvtr(filePrefix + "nvtr.txt");
	//формируем файл rz.txt
	fileXY << rz.size() << endl;
	for (i = 0; i < nPoints; i++)
	{
		fileXY << rz[i].r << " " << rz[i].z << endl;
	}

	//формируем файл nvtr.txt
	nKE = KE.size();
	fileNvtr << nKE << endl;
	for (i = 0; i < nKE; i++)
	{
		for (j = 0; j < 4; j++)
		{
			fileNvtr << KE[i].uzel[j] << " ";
		}
		fileNvtr << endl;
	}
}

void Mka2D_cylindrical::GeneratePortrait()
{
	vector<vector<int>> list;
	list.resize(2);
	int* listbeg = new int[nPoints];
	int listsize = -1;
	int i, j, k, ielem, ind1, ind2, iaddr;

	for (i = 0; i < nPoints; i++)
	{
		listbeg[i] = -1;
	}
	for (ielem = 0; ielem < nKE; ielem++)
	{
		for (i = 0; i < 4; i++)
		{
			k = KE[ielem].uzel[i];
			for (j = i + 1; j < 4; j++)
			{
				ind1 = k;
				ind2 = KE[ielem].uzel[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				if (ind2 == ind1)continue;
				iaddr = listbeg[ind2];
				if (iaddr == -1)
				{
					listsize++;
					listbeg[ind2] = listsize;
					list[0].resize(listsize + 1);
					list[1].resize(listsize + 1);
					list[0][listsize] = ind1;
					list[1][listsize] = -1;
				}
				else
				{
					while (list[0][iaddr] < ind1 && list[1][iaddr] > -1)
						iaddr = list[1][iaddr];
					if (list[0][iaddr] > ind1)
					{
						listsize++;
						list[0].resize(listsize + 1);
						list[1].resize(listsize + 1);
						list[0][listsize] = list[0][iaddr];
						list[1][listsize] = list[1][iaddr];
						list[0][iaddr] = ind1;
						list[1][iaddr] = listsize;
					}
					else if (list[0][iaddr] < ind1)
					{
						listsize++;
						list[0].resize(listsize + 1);
						list[1].resize(listsize + 1);
						list[1][iaddr] = listsize;
						list[0][listsize] = ind1;
						list[1][listsize] = -1;
					}
				}
			}
		}
	}

	ig = new int[nPoints + 1];
	di = new double[nPoints];
	b = new double[nPoints];
	ig[0] = 0;
	vector<int> jgTemp;
	for (i = 0; i < nPoints; i++)
	{
		di[i] = 0;
		b[i] = 0;
		ig[i + 1] = ig[i];
		iaddr = listbeg[i];
		while (iaddr != -1)
		{
			ig[i + 1] = ig[i + 1] + 1;
			jgTemp.push_back(list[0][iaddr]);
			iaddr = list[1][iaddr];
		}
	}
	jg = new int[ig[nPoints]];
	ggl = new double[ig[nPoints]];
	ggu = new double[ig[nPoints]];
	for (i = 0; i < ig[nPoints]; i++)
	{
		ggu[i] = 0;
		ggl[i] = 0;
		jg[i] = jgTemp[i];
	}

	ofstream ig_jg(filePrefix + "ig_jg.txt");
	for (i = 0; i <= nPoints; i++)
	{
		ig_jg << ig[i] << " ";
	}
	ig_jg << endl;
	for (i = 0; i < ig[nPoints]; i++)
	{
		ig_jg << jg[i] << " ";
	}
}

double Mka2D_cylindrical::analiticSolution(Point_cylindrical source)
{
	return 0;
	//return source.z;
}

double Mka2D_cylindrical::Func(int ielem, int i)
{
	return 0;
}

double Mka2D_cylindrical::BasicFunc1d(int num, double ksi)
{
	switch (num)
	{
	case 0: return 1 - ksi;
	case 1: return ksi;
	default:
		cerr << "Error in Basic Function" << endl;
		system("pause");
		exit(1);
	}
}

double Mka2D_cylindrical::BasicFunc2d(int i, double ksi, double eta)
{
	return BasicFunc1d(i % 2, ksi) * BasicFunc1d(i / 2, eta);
}

double Mka2D_cylindrical::difBasicFunc1d(int num, double ksi)
{
	switch (num)
	{
	case 0: return -1;
	case 1: return 1;
	default:
		cerr << "Error in Basic Function" << endl;
		system("pause");
		exit(1);
	}
}

void Mka2D_cylindrical::CreateLocalMatrixs(int ielem)
{
	double hRlocal, hZlocal;
	int i, j;
	double tKoef = sqrt(3. / 5.);
	double integrPoints[3] = {
		0.5,
		(1 + tKoef) / 2,
		(1 - tKoef) / 2
	};
	double tauKoefs[3] = {
		8. / 9.,
		5. / 9.,
		5. / 9.
	};
	hRlocal = rz[KE[ielem].uzel[1]].r - rz[KE[ielem].uzel[0]].r;
	hZlocal = rz[KE[ielem].uzel[2]].z - rz[KE[ielem].uzel[0]].z;
	for (i = 0; i < 4; i++)
	{
		localB[i] = 0;
		for (j = 0; j < 4; j++)
		{
			localMatrix[i][j] = sreda[KE[ielem].numberField].sigma * (helpG1[i][j] * hZlocal * (1.0 + 2.0 * rz[KE[ielem].uzel[0]].r / hRlocal) +
				helpGOseSim[i][j] * hRlocal * hRlocal / hZlocal +
				helpG2[i][j] * 2.0 * rz[KE[ielem].uzel[0]].r * hRlocal / hZlocal) / 12.0;

			//double g = 0, m = 0;
			//for (size_t k = 0; k < 3; k++)
			//{
			//	for (size_t l = 0; l < 3; l++)
			//	{
			//		double d_ksi_i_dr = difBasicFunc1d(i % 2, integrPoints[k])*BasicFunc1d(i / 2, integrPoints[l]);
			//		double d_ksi_i_dz = BasicFunc1d(i % 2, integrPoints[k])*difBasicFunc1d(i / 2, integrPoints[l]);
			//		double d_ksi_j_dr = difBasicFunc1d(j % 2, integrPoints[k])*BasicFunc1d(j / 2, integrPoints[l]);
			//		double d_ksi_j_dz = BasicFunc1d(j % 2, integrPoints[k])*difBasicFunc1d(j / 2, integrPoints[l]);

			//		g += sreda[KE[ielem].numberField].sigma * tauKoefs[k] * tauKoefs[l] * (rz[KE[ielem].uzel[0]].r + integrPoints[k] * hRlocal)*
			//			(d_ksi_i_dr*d_ksi_j_dr + d_ksi_i_dz*d_ksi_j_dz);

			//		//m+=
			//	}
			//}
			//double oldValue = sreda[KE[ielem].numberField].sigma * (helpG1[i][j] * hZlocal * (1.0 + 2.0 * rz[KE[ielem].uzel[0]].r / hRlocal) +
			//	helpGOseSim[i][j] * hRlocal * hRlocal / hZlocal +
			//	helpG2[i][j] * 2.0 * rz[KE[ielem].uzel[0]].r * hRlocal / hZlocal) / 12.0;
			//g *= hRlocal*hZlocal / 4.;
			//localMatrix[i][j] = g;

			//right part
			//double resultRightPart = 0;
			//for (size_t k = 0; k < 3; k++)
			//{
			//	for (size_t l = 0; l < 3; l++)
			//	{
			//		resultRightPart += BasicFunc1d(i % 3, integrPoints[k]) * BasicFunc1d(i / 3, integrPoints[l]) *
			//			Func(Point(
			//				xy[KE[ielem].uzel[0]].x + beta[1] * integrPoints[k] + beta[0] * integrPoints[l] + beta[4] * integrPoints[k] * integrPoints[l],
			//				xy[KE[ielem].uzel[0]].y + beta[3] * integrPoints[k] + beta[2] * integrPoints[l] + beta[5] * integrPoints[k] * integrPoints[l])) *
			//			tauKoefs[k] * tauKoefs[l] * J;
			//	}
			//}
			//localB[i] += resultRightPart / 4. * signAlfa0;
		}
	}
}

void Mka2D_cylindrical::Addition(int ielem)
{
	int ibeg, iend, ind;
	for (int i = 0; i < 4; i++)
	{
		di[KE[ielem].uzel[i]] += localMatrix[i][i];
		b[KE[ielem].uzel[i]] += localB[i];
		for (int j = 0; j < i; j++)
		{
			for (int k = ig[KE[ielem].uzel[i]]; k < ig[KE[ielem].uzel[i] + 1]; k++)
			{
				if (jg[k] == KE[ielem].uzel[j])
				{
					ggl[k] += localMatrix[i][j];
					break;
				}
			}
		}
	}
}

void Mka2D_cylindrical::PrintPlotMatrix(bool flag_simmeric)
{
	int i, j;
	double** APlot = new double*[nPoints];
	for (i = 0; i < nPoints; i++)
	{
		APlot[i] = new double[nPoints];
		for (j = 0; j < nPoints; j++)
		{
			APlot[i][j] = 0;
		}
	}
	if (flag_simmeric)
		for (i = 0; i < nPoints; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggl[j];
			}
		}
	else
		for (i = 0; i < nPoints; i++)
		{
			APlot[i][i] = di[i];
			for (j = ig[i]; j < ig[i + 1]; j++)
			{
				APlot[i][jg[j]] = ggl[j];
				APlot[jg[j]][i] = ggu[j];
			}
		}

	for (i = 0; i < nPoints; i++)
	{
		for (j = 0; j < nPoints; j++)
		{
			testFile << setw(15) << APlot[i][j];
		}
		testFile << endl;
		delete APlot[i];
	}
	testFile << endl;
	delete APlot;

	for (i = 0; i < nPoints; i++)
	{
		testFile << setw(15) << b[i];
	}
	testFile << endl;
	testFile << endl;
}

//индекс точки в массиве xy
int Mka2D_cylindrical::indexRZ(Point_cylindrical source)
{
	return indexRZ(source.r, source.z);
}

int Mka2D_cylindrical::indexRZ(double r, double z)
{
	for (size_t i = 0; i < nPoints; i++)
	{
		if (r == rz[i].r && z == rz[i].z)
			return i;
	}
	cerr << "IndexXY error";
	system("pause");
	exit(1);
	return -1;
}

void Mka2D_cylindrical::Edge1_sim()
{
	ofstream ku1(filePrefix + "ku1.txt");
	vector<int> masEdge1;
	for (int i = 0; i < R.size(); i++)
	{
		if (ku1_down) masEdge1.push_back(indexRZ(R[i], Z[0]));
		if (ku1_up) masEdge1.push_back(indexRZ(R[i], Z[Z.size() - 1]));
	}
	for (int j = 0; j < Z.size(); j++)
	{
		if (ku1_right) masEdge1.push_back(indexRZ(R[R.size() - 1], Z[j]));
		if (ku1_left) masEdge1.push_back(indexRZ(R[0], Z[j]));
	}
	for (int k = 0; k < masEdge1.size(); k++)
	{
		di[k] = 1;
		b[k] = analiticSolution(rz[k]);
		for (int m = ig[k]; m < ig[k + 1]; m++)
		{
			ggl[m] = 0;
		}
		for (int l = 0; l < nPoints; l++)
		{
			for (int m = ig[l]; m < ig[l + 1]; m++)
			{
				if (k == jg[m])
				{
					ggu[m] = 0;
				}
			}
		}
		ku1 << k << '\t' << b[k] << endl;
	}
}

void Mka2D_cylindrical::Edge1_sim_old()
{
	 ofstream ku1(filePrefix + "ku1.txt");
	 vector<int> masEdge1;
	 for (int i = 0; i < R.size(); i++)
	 {
		 if (ku1_down) masEdge1.push_back(indexRZ(R[i], Z[0]));
		 if (ku1_up) masEdge1.push_back(indexRZ(R[i], Z[Z.size() - 1]));
	 }
	 for (int j = 0; j < Z.size(); j++)
	 {
		 if (ku1_right) masEdge1.push_back(indexRZ(R[R.size() - 1], Z[j]));
		 if (ku1_left) masEdge1.push_back(indexRZ(R[0], Z[j]));
	 }
	 for (int iMas = 0; iMas < masEdge1.size(); iMas++)
	 {
		 int k = masEdge1[iMas];
		 di[k] = 1;
		 for (int m = ig[k]; m < ig[k + 1]; m++)
		 {
			 b[jg[m]] -= ggl[m] * analiticSolution(rz[k]);
			 ggl[m] = 0;
		 }
		 b[k] = analiticSolution(rz[k]);
		 for (int l = 0; l < nPoints; l++)
		 {
			 for (int m = ig[l]; m < ig[l + 1]; m++)
			 {
				 if (k == jg[m])
				 {
					 b[l] -= b[k] * ggl[m];
					 ggl[m] = 0;
				 }
			 }
		 }
		 ku1 << k << '\t' << b[k] << endl;
	 }
 }

void Mka2D_cylindrical::Edge2()
{
	ofstream ku2(filePrefix + "ku2.txt");
	double h = 0.001;
	double tetta[2];
	double dh;
	for (int iy = 0; iy < Z.size() - 1; iy++)
	{
		int koef = -1;
		int intIndexThis = indexRZ(R[0], Z[iy]);
		int intIndexNext = indexRZ(R[0], Z[iy + 1]);
		dh = fabs(rz[intIndexNext].z - rz[intIndexThis].z);
		tetta[0] = koef * (analiticSolution(Point_cylindrical(rz[intIndexThis].r + h, rz[intIndexThis].z)) -
			analiticSolution(Point_cylindrical(rz[intIndexThis].r - h, rz[intIndexThis].z))) / (2 * h);
		tetta[1] = koef * (analiticSolution(Point_cylindrical(rz[intIndexNext].r + h, rz[intIndexNext].z)) -
			analiticSolution(Point_cylindrical(rz[intIndexNext].r - h, rz[intIndexNext].z))) / (2 * h);

		//НУЛЕВЫЕ
		tetta[0] = tetta[1] = 0;
		b[intIndexThis] += (2 * tetta[0] + tetta[1]) * dh / 6.0;
		b[intIndexNext] += (tetta[0] + 2 * tetta[1]) * dh / 6.0;
	}
	for (int ix = 0; ix < R.size() - 1; ix++)
	{
		int koef = 1;
		int intIndexThis = indexRZ(R[ix], Z[Z.size() - 1]);
		int intIndexNext = indexRZ(R[ix + 1], Z[Z.size() - 1]);
		dh = fabs(rz[intIndexNext].r - rz[intIndexThis].r);
		tetta[0] = koef * (analiticSolution(Point_cylindrical(rz[intIndexThis].r, rz[intIndexThis].z + h)) -
			analiticSolution(Point_cylindrical(rz[intIndexThis].r, rz[intIndexThis].z - h))) / (2 * h);
		tetta[1] = koef * (analiticSolution(Point_cylindrical(rz[intIndexNext].r, rz[intIndexNext].z + h)) -
			analiticSolution(Point_cylindrical(rz[intIndexNext].r, rz[intIndexNext].z - h))) / (2 * h);

		//НУЛЕВЫЕ
		tetta[0] = tetta[1] = 0;
		b[intIndexThis] += (2 * tetta[0] + tetta[1]) * dh / 6.0;
		b[intIndexNext] += (tetta[0] + 2 * tetta[1]) * dh / 6.0;
	}
}

void Mka2D_cylindrical::GenerateMatrix()
{
	int ielem, i, j;
	for (ielem = 0; ielem < nKE; ielem++)
	{
		CreateLocalMatrixs(ielem);
		Addition(ielem);
	}
	//		PrintPlotMatrix(true);
	//Edge2();
	int sourceIndex = indexRZ(R[0], Z[Z.size() - 1]);
	b[sourceIndex] = power;
	Edge1_sim_old();
	for (size_t i = 0; i < ig[nPoints]; i++)
	{
		ggu[i] = ggl[i];
	}
}

void Mka2D_cylindrical::genNet1d(double startValue, double endValue, double startH, double koefRazriadki, vector<double>& vect)
{
	double direct = 1;
	if (koefRazriadki < 0)
	{
		double tmp = startValue;
		startValue = endValue;
		endValue = tmp;
		direct = -1;
		koefRazriadki *= -1;
	}
	double xCurrent = startValue + startH * direct;
	double h = startH * koefRazriadki;
	vect.push_back(startValue);
	for (size_t i = 0; xCurrent * direct < endValue * direct; i++)
	{
		vect.push_back(xCurrent);
		xCurrent += h * direct;
		h *= koefRazriadki;
	}
	vect.push_back(endValue);
	if (direct < 0)
	{
		int startPos = 0, endPos = vect.size() - 1;
		double tmp;
		while (startPos < endPos)
		{
			tmp = vect[endPos];
			vect[endPos] = vect[startPos];
			vect[startPos] = tmp;
			startPos++;
			endPos--;
		}
	}
}

void Mka2D_cylindrical::insertInVector(double value, vector<double>& vect)
{
	int posInsert = -1;
	for (int i = 0; i < vect.size(); i++)
	{
		if (vect[i] >= value)
		{
			posInsert = i;
			break;
		}
	}
	if (posInsert != -1 && fabs(vect[posInsert] - value) > 1e-9)
	{
		vect.insert(vect.begin() + posInsert, value);
	}
}

//	Умножение матрицы на вектор
void Mka2D_cylindrical::MultMatrixOnVector(double* in, double* out)
{
	int i, j;
	double* out1;
	out1 = new double[nPoints];
	for (i = 0; i < nPoints; i++)
	{
		out1[i] = di[i] * in[i];
		for (j = ig[i]; j < ig[i + 1]; j++)
		{
			out1[i] += ggl[j] * in[jg[j]];
			out1[jg[j]] += ggu[j] * in[i];
		}
	}
	for (i = 0; i < nPoints; i++)
		out[i] = out1[i];
	delete[] out1;
}

double Mka2D_cylindrical::ScalarMult(double* v1, double* v2)
{
	int i;
	double result;
	result = 0;
	for (i = 0; i < nPoints; i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}


void Mka2D_cylindrical::LOS()
{
	int maxiter = 10000, i;
	double alfa, alfachisl, alfaznam, beta, betachisl, betaznam, checkE, epsMSG = 1e-16, startNeviazka;
	double* r = new double[nPoints];
	double* s = new double[nPoints];
	double* z = new double[nPoints];
	double* p = new double[nPoints];
	q = new double[nPoints];
	double* rout = new double[nPoints];
	for (i = 0; i < nPoints; i++)
	{
		s[i] = rout[i] = r[i] = q[i] = z[i] = p[i] = 0;
	}
	MultMatrixOnVector(q, r);
	for (i = 0; i < nPoints; i++)
	{
		r[i] = b[i] - r[i];
		z[i] = r[i];
	}
	MultMatrixOnVector(z, p);
	checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	//startNeviazka = checkE = ScalarMult(r, r);
	for (int iter = 0; iter < maxiter && checkE >= epsMSG; iter++)
	{
		alfachisl = ScalarMult(p, r);
		alfaznam = ScalarMult(p, p);
		alfa = alfachisl / alfaznam;
		for (i = 0; i < nPoints; i++)
		{
			q[i] = q[i] + alfa * z[i];
			r[i] = r[i] - alfa * p[i];
		}
		MultMatrixOnVector(r, rout);
		betachisl = ScalarMult(p, rout);
		betaznam = ScalarMult(p, p);
		beta = -betachisl / betaznam;
		for (i = 0; i < nPoints; i++)
		{
			z[i] = r[i] + beta * z[i];
			p[i] = rout[i] + beta * p[i];
		}
		checkE = sqrt(ScalarMult(r, r) / ScalarMult(b, b));
	}
		double sumPogr = 0;
		double sumU = 0;
		double correctSolution;
		testFile << setw(10) << "r" << setw(10) << "z" << setw(18) << "correctSolution" << setw(18) << "q" << setw(18) << "diff" << endl;
		for (size_t i = 0; i < nPoints; i++)
		{
			correctSolution = analiticSolution(rz[i]);
			testFile << setw(10) << rz[i].r << setw(10) << rz[i].z << setw(18) << correctSolution << setw(18) << q[i] << setw(18) << q[i] - correctSolution << endl;
			//		testFile << setw(16) <<  setw(20) << q[i] << setw(20) << correctSolution << setw(20) << q[i] - correctSolution << setw(5) << i << endl;
			sumPogr += (q[i] - correctSolution) * (q[i] - correctSolution);
			sumU += correctSolution * correctSolution;
		}
		correctSolution = sqrt(sumPogr / sumU);
		testFile << setw(16) <<  "Отн. погрешность: " << correctSolution << endl;
}

void Mka2D_cylindrical::directSolveStraightTask()
{
	GeneratePortrait();
	GenerateMatrix();
	LOS();
}

int Mka2D_cylindrical::findKE(Point point)
{
	return findKE(Point_cylindrical(point.x, point.z));
}

int Mka2D_cylindrical::findKE(Point_cylindrical point)
{
	for (int i = 0; i < KE.size(); i++)
	{
		if (point.r >= rz[KE[i].uzel[0]].r && point.r <= rz[KE[i].uzel[1]].r &&
			point.z >= rz[KE[i].uzel[0]].z && point.z <= rz[KE[i].uzel[2]].z) 
		{
			return i;
		}
	}
	return -1;
	//stringstream ss;
	//ss << "Not found Ke2D for point r=" << point.r << ", z=" << point.z;
	//throw new exception(ss.str().c_str());
}

double Mka2D_cylindrical::solutionInPoint(int iKe, Point_cylindrical target) {
	double result = 0;
	for (int j = 0; j < 4; j++)
	{
		result += q[KE[iKe].uzel[j]] * (1 - fabs(rz[KE[iKe].uzel[j]].r - target.r) / (rz[KE[iKe].uzel[1]].r - rz[KE[iKe].uzel[0]].r))
			* (1 - fabs(rz[KE[iKe].uzel[j]].z - target.z) / (rz[KE[iKe].uzel[2]].z - rz[KE[iKe].uzel[0]].z));
	}
	return result;
}
