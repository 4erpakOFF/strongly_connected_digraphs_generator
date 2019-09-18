#pragma once
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>
#include <random>
#include <ctime>

// Strongly Connected Digraphs
namespace SCDiGraphs  
{
	using AdjacencyMatrix = std::vector< std::vector<short int> >;
	using GraphAndInduced = std::pair< AdjacencyMatrix, std::vector<std::vector<size_t>> >;

	// Strongly Connected Digraphs Generator
	static class Generator
	{
	public:
		static void print(GraphAndInduced);
		static void print(AdjacencyMatrix);

		// ����������� �������. ��� � ������. 
		// ������ true, ���� ���� ������� ������ ��������� ����������� ���� �� ����� ������ ������� �����.
		// ���������� ����������.
		static bool isStronglyConnected(AdjacencyMatrix matrix);   
		static bool isStronglyConnected(GraphAndInduced pair);

		static AdjacencyMatrix generateSimpleGraph(size_t size, double factorMultiplier);
		static AdjacencyMatrix generateSimpleGraph(size_t size);

		// ���������� ����: 
		// First:  ������ ��������� ����, ������� ����������� �� ��������� ����� ����������� ������ ��������� ���������
		// Second: ����� ������ ��������� ���������, �������������� � ���� �������� � ��� ������
		static GraphAndInduced generateGraphAndInducedSubgraphs(size_t minNumOfSubgr, size_t maxNumOfSubgr, size_t minSize, size_t maxSize);
		static GraphAndInduced generateGraphAndInducedSubgraphs(size_t numOfSubgraphs, size_t minSize, size_t maxSize);
		static GraphAndInduced generateGraphAndInducedSubgraphs(size_t numOfSubgraphs, size_t Size);

	private:
		static void printMatrix(AdjacencyMatrix matrix);

		static void AddEdge(AdjacencyMatrix& matrix, size_t fromV, size_t toV);
		static size_t getNumOfEdges(AdjacencyMatrix matrix);
		static bool findEdge(size_t& startI, size_t& startJ, AdjacencyMatrix matrix);

		// �������������� ���������� ���������� ���������� ������, ��� ����������� � ��������������� ����
		static void randomizeVertecies(size_t& fromV, size_t& toV, size_t currentGraphSize);

		static std::vector<AdjacencyMatrix> generateSubgraphs(size_t minNum, size_t maxNum, size_t minSize, size_t maxSize);
		static std::vector<AdjacencyMatrix> generateSubgraphs(size_t numOfGraphs, size_t minSize, size_t maxSize);
		static std::vector<AdjacencyMatrix> generateSubgraphs(size_t numOfGraphs, size_t size);
		static std::vector<AdjacencyMatrix> generateSubgraphs(size_t numOfGraphs);

		//���������� ������, ������� �����������, ���� ��������� ������� ������� ����� � �������������� �����:
		// "array[0] = 3" ������, ��� ������� "0" ������ ���������� "3".
		static  std::vector<size_t> generateNewLabelsForV(AdjacencyMatrix basicMatrix, AdjacencyMatrix additionalMatrix, 
															size_t basicI, size_t basicJ, 
															size_t addI, size_t addJ);

		// ���������� ��� ����� �� ����, �������� ��� ������������� ������ ����������� (����� ����� ���� ������� �� ����������� ��������).
		static void connectTwoGraphsByEdge(AdjacencyMatrix& resultMatrix, AdjacencyMatrix& additionalMatrix,
										   std::vector<std::vector<size_t>>& EmptyListOfStronglyConnectedInducedSubgraphs);

		static GraphAndInduced connectAllSubgraphs(std::vector<AdjacencyMatrix> vectorOfSubgraphs);


		// ������ ���������:

		void justTest(int num, int min, int max);    // ������ �������� �������
		bool oldIsStronglyConnected(AdjacencyMatrix matrix);
		AdjacencyMatrix oldGeneration(int size);   // ��������� ���������� ������ ����������.
	};
}
