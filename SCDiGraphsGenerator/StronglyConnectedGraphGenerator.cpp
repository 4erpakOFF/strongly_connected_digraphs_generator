#include "StronglyConnectedGraphGenerator.h"

using namespace SCDiGraphs;

void Generator::printMatrix(AdjacencyMatrix matrix)
{
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix[i].size(); j++)
			std::cout << matrix[i][j] << "  ";
		std::cout << std::endl;
	}
}

void Generator::print(GraphAndInduced pair)
{
	if (pair.first.empty())
		return;
	if (pair.second.empty())
		return;
	std::cout << "Result Strongly Connected Graph Adjancy Matrix:" << std::endl;
	printMatrix(pair.first);

	std::cout << "Vertices of Induced Strongly Connected Subgraphs:" << std::endl;
	for (size_t i = 0; i < pair.second.size(); i++)
	{
		for (size_t j = 0; j < pair.second[i].size(); j++)
			std::cout << pair.second[i][j] << " ";
		std::cout << std::endl;
	}

}
void Generator::print(AdjacencyMatrix matrix)
{
	return printMatrix(matrix);
}

bool Generator::isStronglyConnected(AdjacencyMatrix matrix)
{
	const short int out = 0, in = 1;
	std::vector<std::vector<short int>> check;  // 0 - out, 1 - in
	check.resize(matrix.size());
	for (size_t i = 0; i < check.size(); i++)
		check[i].resize(2, 0);

	for(size_t i = 0; i < matrix.size(); i++)
		for(size_t j = 0; j < matrix.size(); j++)
			if (matrix[i][j] == 1)
			{
				check[i][out] = 1;
				check[j][in] = 1;
			}
	for (size_t i = 0; i < check.size(); i++)
		if (check[i][in] == 0 || check[i][out] == 0)
			return false;

	return true;
}
bool Generator::isStronglyConnected(GraphAndInduced pair)
{
	return isStronglyConnected(pair.first);
}

void Generator::AddEdge(AdjacencyMatrix& matrix, size_t fromV, size_t toV)
{
	if ((matrix.size() > fromV) && (matrix[fromV].size() > toV) && toV != fromV)
		matrix[fromV][toV] = 1;
	else
		std::cerr << "Impossible to add  " << fromV << "->" << toV << "  edge" << std::endl;
}

size_t Generator::getNumOfEdges(AdjacencyMatrix matrix)
{
	size_t counter = 0;
	for (size_t i = 0; i < matrix.size(); i++)
		for (size_t j = 0; j < matrix[i].size(); j++)
			if (matrix[i][j] == 1)
				counter++;
	return counter;
}

bool Generator::findEdge(size_t& startI, size_t& startJ, AdjacencyMatrix matrix)		//Начиная с указанного элемента, ищет какую-нибудь дугу
{
	size_t i = startI, j = startJ;
	do
	{
		do
		{
			if (matrix[i][j] == 1)
			{
				startI = i;
				startJ = j;
				return true;
			}
			j = (j + 1) % matrix[i].size();
		} while (j != startJ);
		i = (i + 1) % matrix.size();
	} while (i != startI);
	return false;
}

void Generator::randomizeVertecies(size_t& fromV, size_t& toV, size_t currentGraphSize)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<size_t> distr(0, currentGraphSize - 1);
	fromV = distr(gen);
	toV = distr(gen);
	while (fromV == toV)
		toV = distr(gen);
}

// Недостаток алгоритма в том, что он всегда предполагает существование дуги 0 -> 1.
// Отчасти нивелируется тем, что при объединении с другими графами выбирается случайная дуга.
AdjacencyMatrix Generator::generateSimpleGraph(size_t size, double factorMultiplier)
{
	AdjacencyMatrix matrix;
	if (size == 0)
		return matrix;
	if (fabs(factorMultiplier) < 0.02)
		factorMultiplier = 0.02;

	matrix.resize(size);
	for (size_t i = 0; i < size; i++)
	{
		matrix[i].resize(size);
		matrix[i][i] = 0;
	}
	if (size == 1)
	{
		matrix[0][0] = 0;
		return matrix;
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<size_t> distr(2, size);
	size_t cycleSize = distr(gen);			// Размер цикличного уха. Минимальный == 2
	for (size_t i = 0; i < cycleSize - 1; i++)			//создание цикличного уха
		AddEdge(matrix, i, i + 1);
	AddEdge(matrix, cycleSize - 1, 0);

	//Далее идет добавление остальных ушей

	// Чем выше коэффициент, тем больше возможных итераций при создании новых дуг.
	size_t terminationSparsityFactor = (size_t)(3 * size * fabs(factorMultiplier));
	std::uniform_int_distribution<size_t> terminationDistr(0, terminationSparsityFactor - 1);
	
	// Критерий прекращения генерации новых дуг. Зависит от уже добавленных вершин и от factorMultiplier
	bool terminationCriteria = false;					
	size_t currentGraphSize = cycleSize;

	while (!terminationCriteria)
	{
		size_t fromV = 0, toV = 0;
		randomizeVertecies(fromV, toV, currentGraphSize);

		size_t notAddedV = size - currentGraphSize;
		size_t numOfNewV = 0;
		if (notAddedV > 0)
		{
			std::uniform_int_distribution<size_t> newDistr(1, notAddedV);
			numOfNewV = newDistr(gen);
		}

		for (size_t newV = currentGraphSize; newV < currentGraphSize + numOfNewV; newV++)
		{
			AddEdge(matrix, fromV, newV);
			fromV = newV;
		}
		AddEdge(matrix, fromV, toV);
		currentGraphSize += numOfNewV;
		terminationCriteria = currentGraphSize == size && terminationDistr(gen) == 0;
		//выход из цикла <=> будут добавлены все вершины, и из набора случайных чисел [0..(terminationSparsityFactor-1)] выпадет 0.
	}

	std::cerr << "Strongly connected graph of " << size<< " vertices and " << getNumOfEdges(matrix)<< " edges is generated." << std::endl;
	return matrix;
}
AdjacencyMatrix Generator::generateSimpleGraph(size_t size)
{
	return generateSimpleGraph(size, 1.0);
}

std::vector<AdjacencyMatrix> Generator::generateSubgraphs(size_t minNum, size_t maxNum, size_t minSize, size_t maxSize)
{
	std::vector<AdjacencyMatrix> result;
	if (minSize > maxSize || minSize < 2 || maxSize < 2)
	{
		std::cerr << "Incorrect range of subgraphs sizes" << std::endl;
		return result;
	}
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<size_t> numDistr(minNum, maxNum);
	std::uniform_int_distribution<size_t> sizeDistr(minSize, maxSize);

	size_t numOfGraphs = numDistr(gen);
	for (size_t i = 0; i < numOfGraphs; i++)
	{
		auto randomSize = sizeDistr(gen);
		result.push_back(generateSimpleGraph(randomSize));
	}
	return result;
}
std::vector<AdjacencyMatrix> Generator::generateSubgraphs(size_t numOfGraphs, size_t minSize, size_t maxSize)
{
	return generateSubgraphs(numOfGraphs, numOfGraphs, minSize, maxSize);
}
std::vector<AdjacencyMatrix> Generator::generateSubgraphs(size_t numOfGraphs, size_t size)
{
	return generateSubgraphs(numOfGraphs, numOfGraphs, size, size);
}
std::vector<AdjacencyMatrix> Generator::generateSubgraphs(size_t numOfGraphs)
{
	return generateSubgraphs(numOfGraphs, numOfGraphs, 3, 500);
}

std::vector<size_t> Generator::generateNewLabelsForV(AdjacencyMatrix basicMatrix, AdjacencyMatrix additionalMatrix,
													 size_t basicI, size_t basicJ, size_t addI, size_t addJ)
{
	std::vector<size_t> result;                  // ìàññèâ ïîêàçûââàåò, êàêèå áóäóò íîìåðà âåðøèí âòîðîãî ãðàôà â ðåçóëüòèðóþùåì ãðàôå
	result.resize(additionalMatrix.size());
	for (size_t i = 0; i < result.size(); i++)
	{
		if (i == addI)
			result[i] = basicI;
		else if (i == addJ)
			result[i] = basicJ;
		else if (i < addI && i < addJ)
			result[i] = i + basicMatrix.size();
		else if (i > addI && i > addJ)
			result[i] = i + basicMatrix.size() - 2;
		else
			result[i] = i + basicMatrix.size() - 1;
	}

	return result;
}

//Ïåðâûé àðãóìåíò ÿâëÿåòñÿ êàê êàê äîáàâëÿåìîé ìàòðèöåé, òàê è ðåçóëüòèðóþùåé.
void Generator::connectTwoGraphsByEdge(AdjacencyMatrix& resultMatrix, AdjacencyMatrix& additionalMatrix,
												std::vector<std::vector<size_t>>& InducedSubgraphs)
{
	if (additionalMatrix.empty())
		return;
	if (resultMatrix.empty())
	{
		resultMatrix = additionalMatrix;
		std::vector<size_t> newLabels;
		for (size_t i = 0; i < resultMatrix.size(); i++)
			newLabels.push_back(i);
		InducedSubgraphs.push_back(newLabels);
		return;
	}

	size_t resI = 0, resJ = 0, addI = 0, addJ = 0;
	randomizeVertecies(resI, resJ, resultMatrix.size());
	if (!findEdge(resI, resJ, resultMatrix))
		exit(1);
	randomizeVertecies(addI, addJ, additionalMatrix.size());
	if (!findEdge(addI, addJ, additionalMatrix))
		exit(1);

	//ïðîâåðÿåì ïðîòèâîïîëîæíîå íàïðàâëåíèå, è ïðè íåîáõîäèìîñòè äîáàâëÿåì äóãó â ýòîì íàïðàâëåíèè,
	//÷òîáû âïîñëåäñòâèè ïðè ðàçáèåíèè íà ïîðîæäåííûå ïîäãðàôû, îíè áûëè ñèëüíî ñâÿçíûìè.
	if (resultMatrix[resJ][resI] == 1 && additionalMatrix[addJ][addI] == 0)
		AddEdge(additionalMatrix, addJ, addI);
	else if (resultMatrix[resJ][resI] == 0 && additionalMatrix[addJ][addI] == 1)
		AddEdge(resultMatrix, resJ, resI);

	// ìàññèâ ïîêàçûââàåò, êàêèå áóäóò íîìåðà âåðøèí âòîðîãî ãðàôà â ðåçóëüòèðóþùåì ãðàôå
	std::vector<size_t> newLabels = generateNewLabelsForV(resultMatrix, additionalMatrix, resI, resJ, addI, addJ);

	size_t newSize = resultMatrix.size() + additionalMatrix.size() - 2;   // -2 - ïîòîìó ÷òî 2 âåðøèíû âòîðîãî ãðàôà ñîåäèíÿòñÿ ñ âåðøèíàìè ïåðâîãî
	resultMatrix.resize(newSize);
	for (size_t i = 0; i < resultMatrix.size(); i++)
		resultMatrix[i].resize(resultMatrix.size());

	for (size_t i = 0; i < additionalMatrix.size(); i++)
		for (size_t j = 0; j < additionalMatrix[i].size(); j++)
			resultMatrix[newLabels[i]][newLabels[j]] = additionalMatrix[i][j];

	std::sort(newLabels.begin(), newLabels.end());
	InducedSubgraphs.push_back(newLabels);
}

GraphAndInduced Generator::connectAllSubgraphs(std::vector<AdjacencyMatrix> vectorOfSubgraphs)
{
	AdjacencyMatrix resultMatrix;
	std::vector<std::vector<size_t>> inducedSubgraphs;
	inducedSubgraphs.reserve(vectorOfSubgraphs.size());

	if (vectorOfSubgraphs.size() < 1)
	{
		std::cerr << "Incorrect size of subgraphs range" << std::endl;
	}
	else
	{
		size_t capacity = 0;
		for (size_t i = 0; i < vectorOfSubgraphs.size(); i++)
			capacity += vectorOfSubgraphs[i].size() - 2;
		resultMatrix.reserve(capacity);								//âûäåëèëè íóæíûé ðàçìåð ïîä ðåçóëüòèðóþùóþ ìàòðèöó

		for (size_t i = 0; i < vectorOfSubgraphs.size(); i++)
		{
			std::cerr << "\rConnecting: " << i + 1 << "/" << vectorOfSubgraphs.size();
			connectTwoGraphsByEdge(resultMatrix, vectorOfSubgraphs[i], inducedSubgraphs);
		}
		std::cerr << std::endl;
	}
	GraphAndInduced result(resultMatrix, inducedSubgraphs);
	return result;
}

GraphAndInduced Generator::generateGraphAndInducedSubgraphs(size_t minNumOfSubgraphs, size_t maxNumOfSubgraphs, 
															size_t minSize, size_t maxSize)
{
	//èñêëþ÷àåì ÷àñòü îøèáîê ââîäà íåâåðíîãî äèàïàçîíà
	size_t minNum = minNumOfSubgraphs < maxNumOfSubgraphs ? minNumOfSubgraphs : maxNumOfSubgraphs;
	size_t maxNum = maxNumOfSubgraphs > minNumOfSubgraphs ? maxNumOfSubgraphs : minNumOfSubgraphs;
	size_t min = minSize < maxSize ? minSize : maxSize;
	size_t max = maxSize > minSize ? maxSize : minSize;
	std::vector<AdjacencyMatrix> subgraphs = generateSubgraphs(minNum, maxNum, min, max);
	auto result = connectAllSubgraphs(subgraphs);
	return result;
}
GraphAndInduced Generator::generateGraphAndInducedSubgraphs(size_t numOfSubgraphs, size_t minSize, size_t maxSize)
{
	return generateGraphAndInducedSubgraphs(numOfSubgraphs, numOfSubgraphs, minSize, maxSize);
}
GraphAndInduced Generator::generateGraphAndInducedSubgraphs(size_t numOfSubgraphs, size_t size)
{
	return generateGraphAndInducedSubgraphs(numOfSubgraphs, numOfSubgraphs, size, size);
}

//-----------------------------------------------------------------------------------------------------------------
// Ñòàðûå íàðàáîòêè:

bool Generator::oldIsStronglyConnected(AdjacencyMatrix matrix)
{
	clock_t startTime = clock();
	for (size_t i = 0; i < matrix.size(); i++)
	{
		bool isEdgeFind = false;
		for (size_t j = 0; j < matrix[i].size(); j++)
		{
			if (matrix[i][j] == 1)
			{
				isEdgeFind = true;
				break;
			}
		}
		if (isEdgeFind == false) return false;
		for (size_t j = 0; j < matrix[i].size(); j++)
		{
			if (matrix[j][i] == 1)
			{
				isEdgeFind = true;
				break;
			}
		}
		if (isEdgeFind == false) return false;
	}

	clock_t endTime = clock();
	std::cout << "old time: " << ((double)endTime - startTime) / (double)CLOCKS_PER_SEC << std::endl;
	return true;
}
AdjacencyMatrix Generator::oldGeneration(int size)
{
	//srand((unsigned int)time(0));
	AdjacencyMatrix graph;
	graph.resize(size);
	for (int i = 0; i < size; i++)
	{
		graph[i].resize(size);
		graph[i][i] = 0;
	}
	int count = 0;
	for (int i = 0; i < size; i++)
	{
		int j = 0;
		while (i == (j = rand() % size));
			//do nothing
		std::cout << "lol" << std::endl;
		if (graph[i][j] != 1) count++;
		graph[i][j] = 1;
		while (i == (j = rand() % size));
			//do nothing
		std::cout << "lol" << std::endl;
		if (graph[j][i] != 1) count++;
		graph[j][i] = 1;
	}
	std::cout << count << std::endl;
	return graph;
}
void Generator::justTest(int num, int min, int max)
{
	int i = 0;
	while (true)
	{
		std::vector<AdjacencyMatrix> a = generateSubgraphs(num, min, max);
		for (size_t i = 0; i < a.size(); i++) { printMatrix(a[i]); std::cout << std::endl; }

		auto result = connectAllSubgraphs(a);
		if (!isStronglyConnected(result.first)){printMatrix(result.first);break;}
		if (i == 0) break;
		i++;
	}
}
