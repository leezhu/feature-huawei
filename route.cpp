#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <string.h>
#include "route.h"
#include "lib_record.h"
#include "sys/timeb.h"

using namespace std;

/******************************franz' code************************/
#define MAX_NodeSize (660)
#define MAX_EnodeSize (60)
int EdgSizeofInputGragh;
int NodeSizeofInputGragh;
int ENodeSize;
int source;
int destination;
int NodeSizeofSCC;
bool lastisuceedtofindpath;
int isdfsAllENSCC;
bool vs_fstsccdfs;
bool isDfstodest;
bool debug_runingsearchfunction=true;
bool dfsistofindpath;
int preNodesize;
int ResultValueMinCost;
int setMAXValue=11111111;

int RealExistNode[MAX_NodeSize];
int index_node[MAX_NodeSize];
int node_index[MAX_NodeSize];
int ENode[MAX_EnodeSize];
bool isENode[MAX_NodeSize];
int InputGraphMatix[MAX_NodeSize][MAX_NodeSize];
int FlagofInputGraphEdgeMatix[MAX_NodeSize][MAX_NodeSize];
bool SCC_haved[MAX_NodeSize];
int ordnode_ordSCC[MAX_NodeSize];
int ENnumofordSCC[MAX_NodeSize];
bool havedofTopoSort[MAX_NodeSize];
int ordscc_ordtopo[MAX_NodeSize];
int indegofTopoNode[MAX_NodeSize];
int ordtopo_ordscc[MAX_NodeSize];
int ordscchaveNE_seq[MAX_NodeSize];
int seqordscchaveNE_ordtopo[MAX_NodeSize];
bool scc_istoanode4DestnothaveEN[MAX_NodeSize];
bool istoanode4DestnothaveEN[MAX_NodeSize];
int Enodeoutdegree[MAX_NodeSize];
int Enodeindegree[MAX_NodeSize];
int d[MAX_NodeSize];
bool used[MAX_NodeSize];
int firstvs[MAX_NodeSize];
int PreNode[MAX_EnodeSize];
int preNode_seq[MAX_NodeSize];

vector<int>SCC_GraphList[MAX_NodeSize];
vector<int>SCC_rGraphList[MAX_NodeSize];
vector<int>SCC_vsnode;
vector<int>TopoSortEdge[MAX_NodeSize];
vector<int>rTopoSortEdge[MAX_NodeSize];
vector<int>resultpath;

struct timeb StartTime;
struct timeb OverTime;


typedef pair<int,int> P;

struct state
{
    int hop;
    int id;
    int costsum;
    int vsnode[MAX_NodeSize];
    friend bool operator<(const struct state a,const struct state b)
    {
        if(a.costsum<b.costsum)
        {
            return true;
        }
        else
        {
            if(a.costsum>b.costsum)
            {
                return false;
            }
        }
        return false;
    }
} fq,dp,sp;
priority_queue<state>que;
vector<state>preque[MAX_EnodeSize];


/************************************************************ 我的代码 *********************************************/

/********************** 存储路径模板 *********************/
template<class T>
struct PathArc
{
public:
    int arcNum;         //路径中的边编号
    int headNode;	    //边的头节点
    PathArc<T> *preArc;  // 指向路径的前一条边
};

template<class T>
struct PathList
{
public:
    T pathWeight;
    PathArc<T> *pathTail;
    int meetedNodeNum;
    PathList(void)
    {
        meetedNodeNum = 0;
        pathTail = NULL;
        pathWeight = 0;
    }

    ~PathList(void)
    {	 }

    // 像路径中插入新的节点
    void insert(int &arcNum, T &arcWeight, int &headNode, bool & arcTag)
    {
        PathArc<T> *newArc = new PathArc<T>;
        newArc->arcNum = arcNum;
        newArc->preArc = pathTail;
        newArc->headNode = headNode;
        pathTail = newArc;
        pathWeight += arcWeight;
        meetedNodeNum += arcTag;
    }

    // 判断新插入节点是否已经在路径中出现
    bool isLoop(const int &currentNodeNum)
    {

        if (this->pathTail == NULL)
        {
            return 0;
        }
        for (PathArc<T> *arc = this->pathTail; arc != NULL; arc = arc->preArc)
        {
            if (arc->headNode == currentNodeNum)
            {
                return 1;
            }
        }
        return 0;
    }
};


/********************** 存储迭代节点模板 *********************/
template<class T>
struct KspIterNode
{
public:
    PathList<T> path;		// 路径
    int nodeNum;            // 节点数

    KspIterNode(void)
    {	}

    KspIterNode(const T &weight,const int &nodeNum)
    {
        this->nodeNum = nodeNum;
        path.pathWeight = weight;
    }

    ~KspIterNode(void)
    { 	}

    bool operator<(const KspIterNode<T> &node) const
    {
        return  this->path.pathWeight > node.path.pathWeight;
    }
};

/**********************************************************************/


/*************************** 存储图 出边表 入边表 ***********************/
struct GraphArc
{
public:
    int weight;                //图的边的权重
    int tailNode;				// 边的尾节点
    int next;					// 指向相同出发节点的 下一条边
    bool tag;
};


struct InGraphArc
{
public:
    int headNode;
    int next;
    //int weight;
};

struct InGraph
{
public:
    int *node;
    InGraphArc *inArc;
};


struct Graph
{
public:
    int sourse;						//起始节点
    int destination;				//终点
    vector<int> includingSet;		// 必须经过的节点
//	int nodeNum;                    //图的节点数量
    int edgeNum;					//图的边的数量
    int *node;						//存储图的节点数组
    GraphArc *arc;					//存储图的边的数组

    // 字符串转化为数字
    static void strToNum(int &num, char *str,const int &pos0, const int &pos1);
    // 读入图
    static void readGraph(char *graph[5000], const int &edge_num, Graph &linkGraph, InGraph &inLinkGraph);
    // 读入 condition
    static void readCondition(char *condition, Graph &linkGraph);

    // K短路算法
    static void kShortestPath (vector<int> &result, const Graph &linkGraph, const int & terminateTime);
    // 预处理图
    static void preprocess(Graph &linkGraph, InGraph &inLinkGraph, const int & reduceWeight);

};

// 字符串转化为数字
void Graph::strToNum(int &num, char *str, const int &pos0, const int &pos1)
{
    char *tmp = new char[4];
    int i = 0;
    int j = pos0 + 1;
    for(; j < pos1; )
        tmp[i++]=str[j++];
    tmp[i] = 0;

    num = atoi(tmp);
//	cout << "num = " << num << endl;
    delete[] tmp;
}



// 读入图,分别用node存储节点,arc 存储边,linkGraph 表示入边表,inLinkGraph 表示出边表,都是基于头插法建立图
void Graph::readGraph(char *graph[5000], const int &edge_num, Graph &linkGraph, InGraph &inLinkGraph)
{
    linkGraph.edgeNum = edge_num;
    linkGraph.arc = new GraphArc[edge_num];
    inLinkGraph.inArc = new InGraphArc[edge_num];

    linkGraph.node = new int[600];
    inLinkGraph.node = new int[600];
    for (int i = 0; i < 600; i++)
    {
        linkGraph.node[i] = -1;
        inLinkGraph.node[i] = -1;
    }

    int *pos = new int[3];
    int edgeNum = linkGraph.edgeNum - 1; // 最大边编号
    for (int i = 0; i < edgeNum; i++)
    {
        int j = 0;
        int tag = 0;
        while('\n' != graph[i][j])
        {
            if (',' == graph[i][j])
            {
                pos[tag++] = j;
            }
            j++;
        }
        pos[tag] = j;
        int head, tail, weight;
        Graph::strToNum(head, graph[i], pos[0], pos[1]);
        Graph::strToNum(tail, graph[i], pos[1], pos[2]);
        Graph::strToNum(weight, graph[i], pos[2], pos[3]);
        linkGraph.arc[i].weight = weight;
        linkGraph.arc[i].tag = 0;
        //inLinkGraph.inArc[i].weight = weight;
        inLinkGraph.inArc[i].headNode = head;
        linkGraph.arc[i].tailNode = tail;

        linkGraph.arc[i].next = linkGraph.node[head];
        inLinkGraph.inArc[i].next = inLinkGraph.node[tail];
        linkGraph.node[head] = i;
        inLinkGraph.node[tail] = i;
    }

    //  topo  最后一行 没有 回车, 所以单独读入
    int j = 0;
    int tag = 0;
    while(0 != graph[edgeNum][j])
    {
        if (',' == graph[edgeNum][j])
        {
            pos[tag++] = j;
        }
        j++;
    }
    pos[tag] = j;
    int head, tail , weight;
    Graph::strToNum(head, graph[edgeNum], pos[0], pos[1]);
    Graph::strToNum(tail, graph[edgeNum], pos[1], pos[2]);
    Graph::strToNum(weight, graph[edgeNum], pos[2], pos[3]);
    inLinkGraph.inArc[edgeNum].headNode = head;
    linkGraph.arc[edgeNum].tailNode = tail;
    linkGraph.arc[edgeNum].weight = weight;
    linkGraph.arc[edgeNum].tag = 0;
    //inLinkGraph.inArc[edgeNum].weight = weight;
    linkGraph.arc[edgeNum].next = linkGraph.node[head];
    inLinkGraph.inArc[edgeNum].next = inLinkGraph.node[tail];
    linkGraph.node[head] = edgeNum;
    inLinkGraph.node[tail] = edgeNum;

    delete[] pos;
}


// 读入必过节点 condition
void Graph::readCondition(char *condition, Graph &linkGraph)
{
    int i = 0;
    int *pos = new int[2];
    int tag = 0;
    while(tag < 2)
    {

        if (',' == condition[i])
        {
            pos[tag++] = i;
        }
        i++;
    }
    //	cout << "pos[0] = " << pos[0] << " " << " pos[1] = " << pos[1] << endl;
    Graph::strToNum(linkGraph.sourse, condition, -1, pos[0]);
    Graph::strToNum(linkGraph.destination, condition, pos[0], pos[1]);
    delete[] pos;

    int *pos1 = new int[50];

    int tag1 = 0;
    pos1[tag1++] = i - 1;
    while(0 != condition[i])
    {
        if ('|' == condition[i])
        {
            pos1[tag1++] = i;
        }
        i++;
    }
    pos1[tag1] = i;

    int includeNode;
    for (int i = 0; i < tag1; i++)
    {
        //cout << "pos1[i] = " << pos1[i] << endl;
        Graph::strToNum(includeNode, condition, pos1[i], pos1[i + 1]);
        //cout << "includeNode = " << includeNode << endl;
        linkGraph.includingSet.push_back(includeNode);
    }
    delete[] pos1;
}

// 对图做预处理
void Graph::preprocess(Graph & linkGraph, InGraph &inLinkGraph, const int & reduceWeight)
{
    int i = inLinkGraph.node[linkGraph.destination];

    while(i > -1)
    {
        linkGraph.arc[i].tag = 1;
        linkGraph.arc[i].weight = linkGraph.arc[i].weight - reduceWeight;
        i = inLinkGraph.inArc[i].next;
    }

    for(int i = 0; i < linkGraph.includingSet.size(); i++)
    {

        int j = linkGraph.node[linkGraph.includingSet[i]];

        while(j > -1)
        {
            linkGraph.arc[j].weight = linkGraph.arc[j].weight - reduceWeight;
            j = linkGraph.arc[j].next;
        }

        j = inLinkGraph.node[linkGraph.includingSet[i]];
        while(j > -1)
        {
            linkGraph.arc[j].tag = 1;
            linkGraph.arc[j].weight = linkGraph.arc[j].weight - reduceWeight;
            j = inLinkGraph.inArc[j].next;
        }
    }
    delete[] inLinkGraph.node;
    delete[] inLinkGraph.inArc;
}


// K短路算法
void Graph::kShortestPath (vector<int> &result , const Graph &linkGraph, const int & terminateTime)
{
    priority_queue< KspIterNode<int> > ahpHeap;

    KspIterNode<int> currentNode(0, linkGraph.sourse);
    ahpHeap.push(currentNode);

    int minCost = 999;
    PathArc<int> *minPath = NULL;
    int maxMeetedNodeNum = linkGraph.includingSet.size() + 1;
    while (!ahpHeap.empty() )
    {

        currentNode = ahpHeap.top();
        ahpHeap.pop();
        //cout << "currentNode.path.meetedNodeNum = " << currentNode.path.meetedNodeNum << endl;
        if (currentNode.nodeNum == linkGraph.destination)
        {
            struct timeb rawtime;
            ftime(&rawtime);

            static int ms = rawtime.millitm;
            static unsigned long s = rawtime.time;
            int out_ms = rawtime.millitm - ms;
            unsigned long out_s = rawtime.time - s;

            if (out_ms < 0)
            {
                out_s -= 1;
            }
            if(out_s > terminateTime)
            {
                break;
            }

            if(currentNode.path.meetedNodeNum == maxMeetedNodeNum)
            {
                if(minCost > currentNode.path.pathWeight)
                {
                    delete minPath;
                    minCost = currentNode.path.pathWeight;
                    minPath = currentNode.path.pathTail;
                }
                else
                {
                    delete currentNode.path.pathTail;
                }
            }
            continue;
        }

        for (int i = linkGraph.node[currentNode.nodeNum]; i > -1; i = linkGraph.arc[i].next)
        {
            if (currentNode.path.isLoop(linkGraph.arc[i].tailNode))
            {
                continue;
            }
            KspIterNode<int> nextNode;
            nextNode.nodeNum = linkGraph.arc[i].tailNode;
            nextNode.path = currentNode.path;
            nextNode.path.insert(i, linkGraph.arc[i].weight, currentNode.nodeNum, linkGraph.arc[i].tag);
            ahpHeap.push(nextNode);
        }
    }

    if(minPath != NULL)
    {
        for(PathArc<int> *pathArc = minPath; pathArc != NULL; pathArc = pathArc->preArc)
        {
            result.push_back(pathArc->arcNum);
        }
    }
}

/********************************************************************************************************************/

/**********************************franz's code**************************/

void add_edge(int from,int to)
{
    SCC_GraphList[from].push_back(to);
    SCC_rGraphList[to].push_back(from);
}
void DealwithInputGraph()
{
    bool boolindegree;
    bool booloutdegree;
    booloutdegree=false;
    boolindegree=false;
    for(int k=2; k<NodeSizeofInputGragh; k++)
    {
        for(int i=0; i<NodeSizeofInputGragh; i++)
        {
            if(0<InputGraphMatix[k][i])
                booloutdegree=true;
        }
        for(int j=0; j<NodeSizeofInputGragh; j++)
        {
            if(0<InputGraphMatix[j][k])
                boolindegree=true;
        }
        if((false==boolindegree)||(false==booloutdegree))
        {
            RealExistNode[k]=0;
        }
    }
    return ;
}
bool InitJudgeGraphisCorrect()
{
    bool boolindegree;
    bool booloutdegree;
    booloutdegree=false;
    boolindegree=false;
    for(int k=0; k<ENodeSize; k++)
    {
        if(0==RealExistNode[ENode[k]])
        {
            return false;
        }
    }
    for(int j=0; j<NodeSizeofInputGragh; j++)
    {
        if(0<InputGraphMatix[source][j])
            booloutdegree=true;
    }
    for(int i=0; i<NodeSizeofInputGragh; i++)
    {
        if(0<InputGraphMatix[i][destination])
            boolindegree=true;
    }
    if(false==booloutdegree)
        return false;
    if(false==boolindegree)
        return false;
    return true;
}
void GraphArraytoList()
{
    for(int i=0; i<NodeSizeofInputGragh; i++)
    {
        if(0<RealExistNode[i])
        {
            for(int j=0; j<NodeSizeofInputGragh; j++)
            {
                if(j!=i)
                    if(0<RealExistNode[j])
                    {
                        if(0<InputGraphMatix[i][j])
                        {
                            add_edge(i,j);
                        }
                    }
            }
        }
    }
}
void SCCFowardDfs(int v)
{
    if(true==vs_fstsccdfs)
    {
        if(v==destination)
        {
            isDfstodest=true;
        }
        havedofTopoSort[v]=true;
        if(true==isENode[v])
            isdfsAllENSCC++;

    }
    SCC_haved[v]=true;
    for(int i=0; i<SCC_GraphList[v].size(); i++)
    {
        if(!SCC_haved[SCC_GraphList[v][i]])
            SCCFowardDfs(SCC_GraphList[v][i]);
    }
    SCC_vsnode.push_back(v);
}
void SCCreverseDfs(int v,int k)
{
    ordnode_ordSCC[v]=k;
    SCC_haved[v]=true;
    for(int i=0; i<SCC_rGraphList[v].size(); i++)
    {
        if(!SCC_haved[SCC_rGraphList[v][i]]) SCCreverseDfs(SCC_rGraphList[v][i],k);
    }
}
int SCCMainFunction()
{
    isDfstodest=false;
    isdfsAllENSCC=0;
    vs_fstsccdfs=true;
    for(int v=0; v<NodeSizeofInputGragh; v++)
    {
        if(!SCC_haved[v])SCCFowardDfs(v);
        vs_fstsccdfs=false;
    }
    memset(SCC_haved,0,sizeof(SCC_haved));
    int countnum=0;
    for(int i=SCC_vsnode.size()-1; i>=0; i--)
    {
        if(!SCC_haved[SCC_vsnode[i]])SCCreverseDfs(SCC_vsnode[i],countnum++);
    }
    for(int i=0; i<ENodeSize; i++)
    {
        ENnumofordSCC[ordnode_ordSCC[ENode[i]]]++;
    }
    return countnum;
}
bool TopoOrder(int topon)
{
    int top=-1;
    int i;
    int sorti=0;
    for(i=(topon-1); i>=0; i--)
    {
        if(0==indegofTopoNode[i])
        {
            indegofTopoNode[i]=top;
            top=i;
        }
    }
    for(i=0; i<topon; ++i)
    {
        if(top==-1)
        {
            return false;
        }
        int j=top;
        ordtopo_ordscc[sorti]=j;
        sorti++;
        top=indegofTopoNode[top];
        int topofirstenode=false;
        int exmid=0;
        for(int k=0; k<TopoSortEdge[j].size(); ++k)
        {
            if((0==(--indegofTopoNode[TopoSortEdge[j][k]])))
            {
                if(false==topofirstenode)
                {
                    indegofTopoNode[TopoSortEdge[j][k]]=top;
                    top=TopoSortEdge[j][k];
                    topofirstenode=true;
                }
                else
                {
                    if(ENnumofordSCC[top]>0)
                    {
                        exmid=indegofTopoNode[top];
                        indegofTopoNode[TopoSortEdge[j][k]]=exmid;
                        indegofTopoNode[top]=TopoSortEdge[j][k];
                    }
                    else
                    {
                        indegofTopoNode[TopoSortEdge[j][k]]=top;
                        top=TopoSortEdge[j][k];
                    }
                }
            }
        }
    }
}
void getTopoSort4SCC()
{
    int toponodesize=NodeSizeofSCC;
    bool vs_matrixtopograph[MAX_NodeSize][MAX_NodeSize];
    for(int i=0; i<NodeSizeofInputGragh; i++)
    {
        if(true==havedofTopoSort[i])
        {
            for(int j=0; j<NodeSizeofInputGragh; j++)
            {
                if(j!=i)
                    if(true==havedofTopoSort[j])
                    {
                        if(InputGraphMatix[i][j]>0)
                        {
                            if((ordnode_ordSCC[i]!=ordnode_ordSCC[j])&&(false==vs_matrixtopograph[ordnode_ordSCC[i]][ordnode_ordSCC[j]]))//&&(topoedge[ithnode_ithSCC[i]][ithnode_ithSCC[j]]==false))
                            {
                                TopoSortEdge[ordnode_ordSCC[i]].push_back(ordnode_ordSCC[j]);
                                rTopoSortEdge[ordnode_ordSCC[j]].push_back(ordnode_ordSCC[i]);
                                vs_matrixtopograph[ordnode_ordSCC[i]][ordnode_ordSCC[j]]=true;
                                indegofTopoNode[ordnode_ordSCC[j]]++;
                            }
                        }
                    }
            }
        }
    }
    TopoOrder(toponodesize);
    for(int i=0; i<toponodesize; i++)
    {
        ordscc_ordtopo[ordtopo_ordscc[i]]=i;
    }
    seqordscchaveNE_ordtopo[0]=ordscc_ordtopo[source];
    int sequence=1;
    for(int i=0; i<toponodesize; i++)
    {
        if(0<ENnumofordSCC[ordtopo_ordscc[i]])
        {
            ordscchaveNE_seq[ordtopo_ordscc[i]]=sequence;
            seqordscchaveNE_ordtopo[sequence]=i;
            sequence++;
        }
    }
    ordscchaveNE_seq[ordnode_ordSCC[destination]]=sequence;
    seqordscchaveNE_ordtopo[sequence]=ordscc_ordtopo[ordnode_ordSCC[destination]];
    sequence++;
}
void fromdestdfswithoutscchaveEN(int rdv)
{
    scc_istoanode4DestnothaveEN[rdv]=true;
    for(int i=0; i<rTopoSortEdge[rdv].size(); i++)
    {
        if(0<ENnumofordSCC[rTopoSortEdge[rdv][i]])
        {
            scc_istoanode4DestnothaveEN[rTopoSortEdge[rdv][i]]=true;
        }
        if(false==scc_istoanode4DestnothaveEN[rTopoSortEdge[rdv][i]])
            if((rTopoSortEdge[rdv][i]!=ordnode_ordSCC[source])&&(0==ENnumofordSCC[rTopoSortEdge[rdv][i]]))
                fromdestdfswithoutscchaveEN(rTopoSortEdge[rdv][i]);
    }
}
void fromdestdfswithoutEN(int rdv)
{
    istoanode4DestnothaveEN[rdv]=true;
    for(int i=0; i<SCC_rGraphList[rdv].size(); i++)
    {
        if(true==isENode[SCC_rGraphList[rdv][i]])
        {
            istoanode4DestnothaveEN[SCC_rGraphList[rdv][i]]=true;
        }
        if(SCC_rGraphList[rdv][i]!=source)
            if((false==istoanode4DestnothaveEN[SCC_rGraphList[rdv][i]])&&(false==isENode[SCC_rGraphList[rdv][i]]))
                fromdestdfswithoutEN(SCC_rGraphList[rdv][i]);
    }
}

void prepushEnode()
{
    preNodesize=0;
    PreNode[preNodesize]=source;
    preNode_seq[source]=preNodesize;
    preNodesize++;
    for(int i=0; i<ENodeSize; i++)
    {
        PreNode[preNodesize]=ENode[i];
        preNode_seq[ENode[i]]=preNodesize;
        preNodesize++;
    }
    for(int i=0; i<preNodesize; i++)
    {
        fill(d,d+NodeSizeofInputGragh,setMAXValue);
        fill(firstvs,firstvs+NodeSizeofInputGragh,-1);
        priority_queue<P,vector<P>,greater<P> >midque;
        d[PreNode[i]]=0;
        midque.push(P(0,PreNode[i]));
        while(!midque.empty())
        {
            P p=midque.top();
            midque.pop();
            int v=p.second;
            if((v==destination)||((v!=PreNode[i])&&(true==isENode[v]))||(d[v]<p.first))
                continue;
            int u;
            for(int j=0; j<SCC_GraphList[v].size(); j++)
            {
                u=SCC_GraphList[v][j];
                if(d[u]>(d[v]+InputGraphMatix[v][u]))
                {
                    d[u]=d[v]+InputGraphMatix[v][u];

                    firstvs[u]=v;
                    midque.push(P(d[u],u));
                }
            }
        }
        int pre,next,costsum,hop;
        int icci=ordnode_ordSCC[PreNode[i]];
        int iccj;
        int j;
        for(int l=0; l<preNodesize; l++)
        {
            if(i==0&&l==0)
            {
                continue;
            }
            if(l==0)
            {
                j=destination;
            }
            else
                j=PreNode[l];

            if((firstvs[j]!=-1))//&&((true==isENode[j])||j==destination))
            {
                iccj=ordnode_ordSCC[j];
                if((icci==iccj)||((ordscchaveNE_seq[icci]+1)==ordscchaveNE_seq[iccj]))
                {
                    Enodeoutdegree[PreNode[i]]++;
                    Enodeindegree[j]++;
                    pre=j;
                    hop=0;
                    costsum=0;
                    memset(fq.vsnode,-1,sizeof(fq.vsnode));
                    //memset(fq.bit,0,sizeof(fq.bit));
                    next=firstvs[pre];
                    fq.vsnode[pre]=next;
                    //b1=pre/8;
                    //  b2=pre%8;
                    //fq.bit[b1]|=(1<<b2);
                    while(-1!=next)
                    {
                        costsum+=InputGraphMatix[next][pre];
                        pre=next;
                        hop++;
                        next=firstvs[pre];
                        fq.vsnode[pre]=next;
                        /*
                         if(next!=-1){
                             b1=pre/8;
                             b2=pre%8;
                             fq.bit[b1]|=(1<<b2);
                         }
                         */
                    }
                    fq.id=j;
                    fq.hop=hop;
                    fq.costsum=costsum;
                    preque[i].push_back(fq);
                }
            }
        }
    }
}
bool controltime;//=true;
int ansvs[MAX_NodeSize];

void dfsfor14or15(int start,int enodenum,int costsum,int depth)//,unsigned int bit[19])
{
    if(false==controltime)
    {
        return ;
    }
    ftime(&OverTime);
    if(9<=(OverTime.time-StartTime.time))
    {
        controltime=false;
        return ;
    }
    bool isover;
    int seq_preque;
    int franzstart;
    int franzseqstart;
    int franzcost;
    int back_vs[MAX_NodeSize];
    int prequesize;
    int acctop;

    //int mid;
    memcpy(back_vs,ansvs,NodeSizeofInputGragh*4);



    seq_preque=preNode_seq[start];
    prequesize=preque[seq_preque].size();
    //prequesize=(prequesize/10)+1;
    //prequesize=(prequesize/depth)+1;
    preNodesize=(preNodesize*2/3)+1;
    for(int i=0; i<prequesize; i++)
    {
        if(true==controltime)
        {
            isover=true;
            dp=preque[seq_preque][i];
            if(dp.id!=destination)
            {
                for(int j=0; j<NodeSizeofInputGragh; j++)
                {
                    if((-1!=ansvs[j])&&(-1!=(dp.vsnode[j])))
                    {
                        isover=false;
                        break;
                    }
                    //ansvs[j]=dp.vsnode[j];
                }
                if(isover)
                {
                    if((enodenum+1)==ENodeSize)
                    {
                        if(istoanode4DestnothaveEN[dp.id]&&scc_istoanode4DestnothaveEN[ordnode_ordSCC[dp.id]])
                        {
                            franzstart=dp.id;
                            franzseqstart=preNode_seq[franzstart];
                            franzcost=costsum+dp.costsum;
                            for(int j=0; j<preque[franzseqstart].size(); j++)
                            {
                                sp=preque[franzseqstart][j];
                                if(sp.id==destination)
                                {
                                    break;
                                }

                            }

                            if(sp.vsnode[destination]!=-1)
                            {
                                acctop=dp.hop+sp.hop;
                                for(int j=0; j<NodeSizeofInputGragh; j++)
                                {
                                    if(-1!=(dp.vsnode[j]))
                                    {
                                        ansvs[j]=dp.vsnode[j];
                                        acctop--;
                                    }
                                    if(-1!=(sp.vsnode[j]))
                                    {
                                        ansvs[j]=sp.vsnode[j];
                                        acctop--;
                                    }
                                    if(0==acctop)
                                        break;
                                }
                                //controltime=false;
                                franzcost+=sp.costsum;
                                if(ResultValueMinCost>franzcost)
                                {
                                    ResultValueMinCost=franzcost;
                                    resultpath.clear();
                                    resultpath.push_back(destination);
                                    for(int pathi=ansvs[destination]; pathi!=-1; pathi=ansvs[pathi])
                                    {
                                        resultpath.push_back(pathi);
                                    }
                                    dfsistofindpath=true;
                                }
                            }
                        }
                    }
                    else
                    {
                        acctop=dp.hop;
                        for(int j=0; j<NodeSizeofInputGragh; j++)
                        {
                            if(-1!=(dp.vsnode[j]))
                            {
                                ansvs[j]=dp.vsnode[j];
                                //b1=j/8;
                                // b2=j%8;
                                //bit[b1]|=(1<<b2);
                                acctop--;
                            }
                            if(0==acctop)
                                break;
                        }

                        dfsfor14or15(dp.id,enodenum+1,costsum+dp.costsum,(depth+2));
                        memcpy(ansvs,back_vs,NodeSizeofInputGragh*4);
                        //memcpy(bit,bt,19*4);
                    }
                }
            }
        }
    }
    return ;
}
void bfsfor14or15()
{
    ResultValueMinCost=setMAXValue;
    prepushEnode();

    for(int i=1; i<preNodesize; i++)
        sort(preque[i].begin(),preque[i].end());

    // unsigned int bit[19];
    fill(ansvs,ansvs+NodeSizeofInputGragh,-1);
    //fill(bit,bit+19,0);
    controltime=true;
    //unsigned int bit1,bit2;
    //bit1=source/8;
    //bit2=source%8;
    //bit[bit1]|=(1<<bit2);
    dfsfor14or15(source,0,0,2);//,bit);
    //cout<<"minimumcost:    "<<minimumcost<<endl;

}

/************************************************************************/
void search_route(char *topo[5000], int edge_num, char *demand)
{
    if (edge_num > 2300)
    {
        ftime(&StartTime);
        int midInput;
        int ithtopoSize;
        unsigned short franzresult[MAX_NodeSize];
        int demandStrLen = strlen(demand);
        int EdgeFlag,inEdgeFlag,outEdgeflag,EdgeWeight;

        /********************************************/
        NodeSizeofInputGragh=0;
        EdgSizeofInputGragh=0;
        ENodeSize=0;
        lastisuceedtofindpath=true;
        memset(index_node,-1,sizeof(index_node));
        memset(node_index,-1,sizeof(node_index));

        int demi=0;
        for(midInput=0; demi<demandStrLen; demi++)
        {
            if(!((demand[demi]>='0')&&(demand[demi]<='9')))
            {
                demi++;
                break;
            }
            midInput=midInput*10+(demand[demi]-'0');
        }
        //source = mInput;
        destination = midInput;
        for(midInput=0; demi<demandStrLen; demi++)
        {
            if(!((demand[demi]>='0')&&(demand[demi]<='9')))
            {
                demi++;
                break;
            }
            midInput=midInput*10+(demand[demi]-'0');
        }
        //destination = mInput;
        source = midInput;

        index_node[source]=NodeSizeofInputGragh;
        node_index[NodeSizeofInputGragh]=source;
        source=NodeSizeofInputGragh;

        RealExistNode[NodeSizeofInputGragh]=1;
        NodeSizeofInputGragh++;

        index_node[destination]=NodeSizeofInputGragh;
        node_index[NodeSizeofInputGragh]=destination;
        destination=NodeSizeofInputGragh;

        RealExistNode[NodeSizeofInputGragh]=1;
        NodeSizeofInputGragh++;

        for(midInput=0; demi<=demandStrLen; demi++)
        {
            if(demand[demi]!='\n')
            {
                if((demand[demi]==',')||(demand[demi]=='|')||(demi==(demandStrLen)))
                {
                    index_node[midInput]=NodeSizeofInputGragh;
                    node_index[NodeSizeofInputGragh]=midInput;
                    ENode[ENodeSize]=NodeSizeofInputGragh;
                    isENode[NodeSizeofInputGragh]=true;
                    RealExistNode[NodeSizeofInputGragh]=1;
                    midInput=0;
                    NodeSizeofInputGragh++;
                    ENodeSize++;
                    if(demi==demandStrLen)
                    {
                        break;
                    }
                }
                else
                {
                    midInput=midInput*10+(demand[demi]-'0');
                }
            }
            if((demand[demi]!='\n')&&(demi==demandStrLen))
                break;
        }

        for(int i=0; i<edge_num; i++)
        {
            int j=0;
            midInput=0;
            ithtopoSize=strlen(topo[i]);
            for(; j<ithtopoSize; j++)
            {
                if(!((topo[i][j]>='0')&&(topo[i][j]<='9')))
                {
                    j++;
                    break;
                }
                midInput=midInput*10+(topo[i][j]-'0');
            }
            EdgeFlag=midInput;

            midInput=0;
            for(; j<ithtopoSize; j++)
            {
                if(!((topo[i][j]>='0')&&(topo[i][j]<='9')))
                {
                    j++;
                    break;
                }
                midInput=midInput*10+(topo[i][j]-'0');
            }
            if(index_node[midInput]==-1)
            {
                index_node[midInput]=NodeSizeofInputGragh;
                node_index[NodeSizeofInputGragh]=midInput;
                RealExistNode[NodeSizeofInputGragh]=1;
                NodeSizeofInputGragh++;
            }
            inEdgeFlag=index_node[midInput];

            midInput=0;
            for(; j<ithtopoSize; j++)
            {
                if(!((topo[i][j]>='0')&&(topo[i][j]<='9')))
                {
                    j++;
                    break;
                }
                midInput=midInput*10+(topo[i][j]-'0');
            }

            if(index_node[midInput]==-1)
            {
                index_node[midInput]=NodeSizeofInputGragh;
                node_index[NodeSizeofInputGragh]=midInput;
                RealExistNode[NodeSizeofInputGragh]=1;
                NodeSizeofInputGragh++;

            }
            outEdgeflag=index_node[midInput];
            midInput=0;
            for(; j<ithtopoSize; j++)
            {
                if(!((topo[i][j]>='0')&&(topo[i][j]<='9')))
                {
                    j++;
                    break;
                }
                midInput=midInput*10+(topo[i][j]-'0');
            }
            EdgeWeight=midInput;
            int exchange=inEdgeFlag;
            inEdgeFlag=outEdgeflag;
            outEdgeflag=exchange;
            if((inEdgeFlag!=destination)&&(outEdgeflag!=source)&&(inEdgeFlag!=outEdgeflag))
            {

                if(0==InputGraphMatix[inEdgeFlag][outEdgeflag])
                {
                    InputGraphMatix[inEdgeFlag][outEdgeflag]=EdgeWeight;
                    FlagofInputGraphEdgeMatix[inEdgeFlag][outEdgeflag]=EdgeFlag;
                    EdgSizeofInputGragh++;

                }
                else
                {
                    if(InputGraphMatix[inEdgeFlag][outEdgeflag]>EdgeWeight)
                    {
                        InputGraphMatix[inEdgeFlag][outEdgeflag]=EdgeWeight;
                        FlagofInputGraphEdgeMatix[inEdgeFlag][outEdgeflag]=EdgeFlag;
                    }
                }

            }
        }

        DealwithInputGraph();
        if(lastisuceedtofindpath==true)
            lastisuceedtofindpath=InitJudgeGraphisCorrect();
        GraphArraytoList();
        NodeSizeofSCC=SCCMainFunction();
        if((isDfstodest==false)||(isdfsAllENSCC!=ENodeSize))
            lastisuceedtofindpath=false;
        getTopoSort4SCC();
        fromdestdfswithoutEN(destination);
        fromdestdfswithoutscchaveEN(ordnode_ordSCC[destination]);
        bfsfor14or15();
        lastisuceedtofindpath=dfsistofindpath;
        if(true==lastisuceedtofindpath)
        {
            int j=0;
            for(int i=resultpath.size()-1; i>=1; i--,j++)
            {
                franzresult[j]=FlagofInputGraphEdgeMatix[resultpath[i]][resultpath[i-1]];
            }
            //for (int i = 0; i < j; i++)
            //record_result(franzresult[i]);
            for (int i = (j-1); i>=0; i--)
                record_result(franzresult[i]);
        }
        /*****************************************************************************************/
    }
    else
    {
        int reduceWeight;    
        int terminateTime;
        Graph linkGraph;
        InGraph inLinkGraph;
        Graph::readGraph(topo, edge_num, linkGraph, inLinkGraph);
        Graph::readCondition(demand, linkGraph);

        if(edge_num <= 100)
        {
            reduceWeight = 0;
            terminateTime = 0;
        }
        if(edge_num > 100 && edge_num <= 300)
        {
            reduceWeight = 10;
            terminateTime = 0;
        }
        if(edge_num > 300 && edge_num <= 430)
        {
            reduceWeight = 8;
            terminateTime = 8;
        }
        if(edge_num > 430 && edge_num <= 580)
        {
            reduceWeight = 10;
            terminateTime = 8;
        }
        if(edge_num > 580 && edge_num <= 1110)
        {
            reduceWeight = 7;
            terminateTime = 1;
        }
        if(edge_num > 1100 && edge_num <= 1200)
        {
            reduceWeight = 10;
            terminateTime = 18;
        }
        if(edge_num > 1200 && edge_num <= 2000 && linkGraph.includingSet.size() <= 23)
        {
            reduceWeight = 6;
            terminateTime = 5;
        }
        if(edge_num > 1200 && edge_num <= 2000  && linkGraph.includingSet.size() > 23)
        {
            reduceWeight = 8;
            terminateTime = 1;
        }
        if(	edge_num > 2000 && edge_num <= 2300)
        {
            reduceWeight = 8;
            terminateTime = 3;
        }

        Graph::preprocess(linkGraph, inLinkGraph, reduceWeight);
        vector<int> result;
        Graph::kShortestPath (result, linkGraph, terminateTime);

        int j = result.size();
        while(j > 0)
        {
            record_result(result[--j]);
        }
    }
}
