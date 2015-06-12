#include "basic/pregel-dev.h"
#include "utils/type.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>

using namespace std;

//input line format: vertexID \t x y
//output line format: vertexID x y centroidID centroidX centroidY

//Global variables
uint u32K;

typedef struct 
{
	double dbx;
	double dby;
} coordinate;

ibinstream & operator<<(ibinstream & m, const coordinate & stCoordnt)
{
	m << stCoordnt.dbx;
	m << stCoordnt.dby;
	return m;
}

obinstream & operator>>(obinstream & m, coordinate & stCoordnt)
{
	m >> stCoordnt.dbx;
	m >> stCoordnt.dby;
	return m;
}
//=============================================

typedef struct 
{
	int iClusterID;
	coordinate stCoordnt;
} KMeansValue;

ibinstream & operator<<(ibinstream & m, const KMeansValue & v)
{
	m << v.iClusterID;
	m << v.stCoordnt;
	return m;
}

obinstream & operator>>(obinstream & m, KMeansValue & v)
{
	m >> v.iClusterID;
	m >> v.stCoordnt;
	return m;
}
//==============================================

class KmeansVertex: public Vertex<VertexID, KMeansValue, VertexID>
{
private:
	double computeDist(coordinate stCentroidCoordnt)
	{
		return sqrt(pow((value().stCoordnt.dbx - stCentroidCoordnt.dbx), 2) + pow((value().stCoordnt.dby - stCentroidCoordnt.dby), 2));		
	}

public:
	virtual void compute(MessageContainer & messages)
	{
		if (step_num() == 1)
		{
			activate();
		}
		else
		{
			int iLastClusterID = value().iClusterID;

			vector<coordinate> *pvctAggCentroids = (vector<coordinate> *)getAgg();
			double dbMinDist = 0;

			double dbCurDist = computeDist((*pvctAggCentroids)[0]);
			dbMinDist = dbCurDist;
			value().iClusterID = 0;

			for (int ii = 1; ii < u32K; ii++)
			{
				dbCurDist = computeDist((*pvctAggCentroids)[ii]);
				if (dbCurDist < dbMinDist)
				{
					dbMinDist = dbCurDist;
					value().iClusterID = ii;
				}
			}

			if (iLastClusterID == value().iClusterID)
			{
				vote_to_halt();
			}
			else
			{
				wakeAll();
			}
		}
	}

};
//=============================================

struct partialSum
{
	int iClusterSize;
	coordinate stCoordntSum;
};

ibinstream & operator<<(ibinstream & m, const partialSum & p)
{
	m << p.iClusterSize;
	m << p.stCoordntSum;
}

obinstream & operator>>(obinstream & m, partialSum & p)
{
	m >> p.iClusterSize;
	m >> p.stCoordntSum;
}

class KmeansAgg: public Aggregator<KmeansVertex, vector<partialSum>, vector<coordinate> >
{
private:
	vector<partialSum> vctPartials;
	vector<coordinate> vctFinals;

	void initializeVectors()
	{
		partialSum stInitSum;
		coordinate stInitCoordnt;
		stInitCoordnt.dbx = 0;
		stInitCoordnt.dby = 0;

		stInitSum.iClusterSize = 0;
		stInitSum.stCoordntSum = stInitCoordnt;

		for (int ii = 0; ii < u32K; ii++)
		{
			vctPartials.push_back(stInitSum);
			vctFinals.push_back(stInitCoordnt);
		}
	}

	void resetPartialSum()
	{
		vctPartials.empty();

		partialSum stInitSum;
		coordinate stInitCoordnt;
		stInitCoordnt.dbx = 0;
		stInitCoordnt.dby = 0;

		stInitSum.iClusterSize = 0;
		stInitSum.stCoordntSum = stInitCoordnt;

		for (int ii = 0; ii < u32K; ii++)
		{
			vctPartials.push_back(stInitSum);
		}
	}


	void calulateCentroids()
	{
		for (int ii = 0; ii < u32K; ii++)
		{
			if (vctPartials[ii].iClusterSize == 0)
			{
				vctFinals[ii].dbx = 0;
				vctFinals[ii].dby = 0;
			}
			else
			{
				vctFinals[ii].dbx = vctPartials[ii].stCoordntSum.dbx / vctPartials[ii].iClusterSize;
				vctFinals[ii].dby = vctPartials[ii].stCoordntSum.dby / vctPartials[ii].iClusterSize;
			}
		}
	}

public:
	virtual void init()
	{
		//Initilize the vectors
		if (step_num() == 1)
		{
			initializeVectors();	
		}
		else
		{
			resetPartialSum();
		}
		

	}

	virtual void stepPartial(KmeansVertex *pv)
	{
		int iClusterID = pv->value().iClusterID;
		vctPartials[iClusterID].iClusterSize += 1;
		vctPartials[iClusterID].stCoordntSum.dbx += pv->value().stCoordnt.dbx;
		vctPartials[iClusterID].stCoordntSum.dby += pv->value().stCoordnt.dby;
	}

	virtual void stepFinal(vector<partialSum> *pvctPartials)
	{
		for (int ii = 0; ii < u32K; ii++)
		{
			vctPartials[ii].iClusterSize += (*pvctPartials)[ii].iClusterSize;
			vctPartials[ii].stCoordntSum.dbx += (*pvctPartials)[ii].stCoordntSum.dbx;
			vctPartials[ii].stCoordntSum.dby += (*pvctPartials)[ii].stCoordntSum.dby;
		}
	}

	virtual vector<partialSum>* finishPartial()
	{
		return &vctPartials;
	}

	virtual vector<coordinate>* finishFinal()
	{

		calulateCentroids();
		
		return &vctFinals;
	}
};
//================================================

class kMeansWorker: public Worker<KmeansVertex, KmeansAgg>
{
	char buf[100];

public:
	virtual KmeansVertex* toVertex(char *pline)
	{
		char *pch;
		pch = strtok(pline, "\t");
		KmeansVertex *pv = new KmeansVertex;
		pv->id = atoi(pch);
		pch = strtok(NULL, " ");
		pv->value().stCoordnt.dbx = atof(pch);
		pch = strtok(NULL, " ");
		pv->value().stCoordnt.dby = atof(pch);

		// Randomly assign the vertex to a cluster
		pv->value().iClusterID = rand() % ::u32K;

		return pv;
	}

	virtual void toline(KmeansVertex *pv, BufferedWriter & writer)
	{
		vector<coordinate> *pvctAggCentroids = (vector<coordinate> *)getAgg();
		int iClusterID = pv->value().iClusterID;
		coordinate stCurrCentroid = (*pvctAggCentroids)[iClusterID];

		// VertexID 	(x, y)	clusterID 	(x, y)	
		sprintf(buf, "%d\t%d %d\t%d\t%d %d\n", pv->id, (int)pv->value().stCoordnt.dbx, (int)pv->value().stCoordnt.dby, iClusterID, (int)stCurrCentroid.dbx, (int)stCurrCentroid.dby);
		writer.write(buf);
	}
};

void pregel_kmeans(string in_path, string out_path, uint k)
{
	::u32K = k;
	srand(time(NULL));
	WorkerParams param;
	param.input_path = in_path;
	param.output_path = out_path;
	param.force_write = true;
	param.native_dispatcher = false;
	kMeansWorker worker;
	KmeansAgg agg;
	worker.setAggregator(&agg);
	worker.run(param);
}
