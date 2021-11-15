package cpcs324_project;
import java.util.Collections;
import javafx.util.Pair;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Scanner;


public class CPCS324_Project {


    ////////////////////////////////class Edge ///////////////////////////////////////////////
    static class Edge implements Comparable<Edge>{
       //declare data field.
        int src;
        int dest;
        int weight;

        // constructor with source, destination, and weight parameter.
        public Edge(int src, int dest, int weight) {
            this.src = src;
            this.dest = dest;
            this.weight = weight;
        }
       
         public int compareTo (Edge E){
       return this.weight-E.weight;
   }
    } 
    ////////////////////////// class HeapNode  ////////////////////////////////////////////
    static class HeapN {

        int v;//declare vertex veraible.
        int k;// declare key verailble.
    }
    ///////////////////////////class result set ////////////////////////////////////////////
    static class Result {
        int parent;//declare parent veraible.
        int weight; //declare weight veraible.
    }   
    /////////////////////////class garph ///////////////////////////////////////////////////////
    // we added the Methods to Graph class.
    static class Graph {
        int vertices;
        int edges;
        LinkedList<Edge>[] adjacencylist; 
       static LinkedList<Edge> e=new LinkedList<Edge>();
        // constructor with vertices, and edges parameter.
        Graph(int vertices, int edges) {
            this.vertices = vertices;
            this.edges = edges;
            adjacencylist = new LinkedList[vertices];
            
            // to initialize adjacency lists.
            for (int k = 0; k < vertices; k++) {
                adjacencylist[k] = new LinkedList<>();
            }
        }
        // Mathod to create Edges in make graph Method        
        public void addEdge(int src, int dest, int weight) {
            Edge edge = new Edge(src, dest, weight);
            adjacencylist[src].addFirst(edge);//add the first edge.
               e.add(edge);
            edge = new Edge(dest, src, weight);
            adjacencylist[dest].addFirst(edge); //undirected graph
             e.add(edge);
        }
        // Prim's Algorithm in pq
        public void primpq() {
            boolean[] minimum_spanning_tree = new boolean[vertices];
            Result[] result = new Result[vertices];
            //keys used to store the key
            int[] KEY = new int[vertices];   
            //Initialize keys and result for all the vertices
            int i=0;
            while(i < vertices){
                KEY[i] = Integer.MAX_VALUE;
                result[i] = new Result();
               i++; //increment by 1 
            }
            
            //create priority queue and we override the comparator
            PriorityQueue<Pair<Integer, Integer>> priorityQ = new PriorityQueue<>(vertices, new Comparator<Pair<Integer, Integer>>() {
               
                
                @Override
                public int compare(Pair<Integer, Integer> P_1, Pair<Integer, Integer> P_2) {
                    int K_1 = P_1.getKey();
                    int K_2 = P_2.getKey();
                    return K_1 - K_2;
                }
            }
            );
            //create the pair 
            KEY[0] = 0;
            Pair<Integer, Integer> P_0 = new Pair<>(KEY[0], 0);
            //add it to  priorityQ
            priorityQ.offer(P_0);
            result[0] = new Result();
            result[0].parent = -1;
                     
            //loop until queue is empty 
            while (priorityQ.isEmpty()==false) {
                Pair<Integer, Integer> pair = priorityQ.poll();
                //extracted vertex
                int vertex = pair.getValue();
                minimum_spanning_tree[vertex] = true;

                //iterate through all the adjacent vertices and update the keys
                LinkedList<Edge> edges = adjacencylist[vertex];
                int j=0;
                while(j< edges.size()){
                    Edge e = edges.get(j);
                    //if edge destination not in minimum_spanning_tree
                    if (!minimum_spanning_tree[e.dest]) {
                        int dest = e.dest;
                        int newK = e.weight; //new key equal edge weight
                        //check if destination >new key 
                        if (KEY[dest] > newK) {
                            //add to pq
                            Pair<Integer, Integer> p = new Pair<>(newK, dest);
                            priorityQ.offer(p);
                            //update the result 
                            result[dest].parent = vertex;
                            result[dest].weight = newK;
                            //update the key[]
                            KEY[dest] = newK;
                        }
                    }
                   j++;//increment by 1 
                }
            }
            //print  minimum_spanning_tree
            printMST(result);
        }
        ///////////////////////////////////////// Prim's Algorithm using Minheap////////////////////////////////////////////////////////   
        public void MinHeap_Based_Prime() {
                
           // keeping track of the vertices which are currently in min heap using insideHeap[]  
            boolean[] insideHeap = new boolean[vertices];
            Result[] FinalSet = new Result[vertices];
            
            
           //Create KEY[] to keep track of key value for each vertex. (keep updating it as heapNode key for each vertex)
            int[] KEY = new int[vertices];
            
            
           // Create a heapN for each vertex which will store two information vertex and key
           
            HeapN[] heapN = new HeapN[vertices];
            for (int j = 0; j < vertices; j++) {
                heapN[j] = new HeapN();
                heapN[j].v = j;
                heapN[j].k = Integer.MAX_VALUE;
                FinalSet[j] = new Result();
                FinalSet[j].parent = -1;  
                insideHeap[j] = true;
                KEY[j] = Integer.MAX_VALUE;
            }

            
            //For each heapNode, Initialize key as MAX_VAL except the heapNode for the first vertex for which key will 0.
            heapN[0].k = 0;

            
            // Create min Heap of size = no of vertices.
            Min_heap min_Heap = new Min_heap(vertices);
            
            
          // Insert all the heapNodes into min heap. inHeap[v] = true for all vertices.
            for (int j = 0; j < vertices; j++) {
                min_Heap.insert(heapN[j]);
            }

            //while min_Heap is not empty
            while (!min_Heap.isEmpty()) {
                
                //Extract the min node from the heap, say it vertex u and add it to the MST.
                HeapN eNode = min_Heap.eMin();

                //extracted vertex
                int eVertex = eNode.v;
                insideHeap[eVertex] = false;

                //Loop for all the adjacent vertices
                LinkedList<Edge> list_E = adjacencylist[eVertex];
                for (int j = 0; j < list_E.size(); j++) {
                    Edge edge = list_E.get(j);
                    
                    //only if edge dest is a present in the heap.
                    if (insideHeap[edge.dest]) {
                        int dest = edge.dest;
                        int newK = edge.weight;
                        
                        //check if the newK is less  than existing key then update 
                        if (KEY[dest] > newK) {
                            dec_Key(min_Heap, newK, dest);
                            
                            //update the parent node.
                            FinalSet[dest].parent = eVertex;
                            FinalSet[dest].weight = newK;
                            KEY[dest] = newK;
                        }
                    }
                }
            }
        
            
            //print MST
            printMST(FinalSet);
        }       
        /////////////////////////////////////////// decrease the key method ///////////////////////////////////////////
        public void dec_Key(Min_heap min_Heap, int newKey, int vertex) {

            //get the index.
            int ind = min_Heap.ind[vertex];

            //get the Node
            HeapN Node = min_Heap.min_H[ind];
            Node.k = newKey;// update the value
            min_Heap.bubbleU(ind);
        }      
        ///////////////////////////////////////print the MST //////////////////////////////////////////////////////////      
        public void printMST(Result[] FinalS) {
            int min_cost = 0;
            System.out.println("Minimum Spanning Tree: ");
            for (int j = 1; j < vertices; j++) {
                min_cost += FinalS[j].weight;
            }
            System.out.println("Total cost: " + min_cost);
        }
       
        ////////////////////////////////////// Kruskal's Algorithm //////////////////////////////////////////////////////////
        
  public static void kruskalMST(  int vertices) {
            
        LinkedList<Edge> mst = new LinkedList<Edge>();//create mst linkedlist 
        LinkedList<Integer> parents = new LinkedList<Integer>();//create parents linkedlist.
        LinkedList<Integer> size = new LinkedList<Integer>();//create size linked list.
        Collections.sort(e);//to sorting linkedlist that storing edges.
        int min=0;//initialize min with 0.
        int ecounter = 0, edgeTaken = 1;//initialize ecounter by 0, and edgeTaking by 1.
        
              //initially the vertices are a set in 
              //themselves.  they are the parent of their
              //own set, Also the size of their set is 1
         
        for (int v =0; v < e.size(); v++) {
           parents.add(v, v);
          size.add(v, 1);
        }
        //loop while will continue until the edgeTaken<=vertices-1.
        while (edgeTaken <= vertices-1) {
            Edge edge = e.get(ecounter);
            ecounter++;//increment by 1.
            if (isCyclic(edge.src, edge.dest, parents)) {//invoke isCyclic method to check the edge will not create cycle
                                                                   //in minmum spanning tree.
                continue;
            }
              /*
              for combining the edge into the minmum spanning tree ,first
             need to find the parents of both the vertices, 
             and then the put combine the smaller set with 
             larger set
             */
            //call union Method.
            union(find(edge.src, parents),
                   //call find Method.
                    find(edge.dest, parents), parents, size);
            mst.add(new Edge(edge.weight, edge.src, edge.dest));
            edgeTaken++;
            min+=edge.weight;//add weight to get mincost.
        }
        
            System.out.println("Minimum Spanning Tree: ");
             System.out.println("Total cost: " + min);    
  }        
        //////////////////////////////////////////find method //////////////////////////////////////////////////////////////
        
      public static int find(int u, LinkedList<Integer> parents) {
          /*if the parent of any vertex is the vertex itself 
          then it is the parent of the the vertex of the
           current edge.
         */
        if (parents.get(u) == u) {
            return u;
        } else {
            //else
            parents.set(u, find(parents.get(u), parents));
            return parents.get(u);
        }
    }

        //////////////////////////////////////////union method //////////////////////////////////////////////////////////////
        
       public static void union(int u, int v, LinkedList<Integer> parents,LinkedList<Integer> size) {
         //find the parent of both the vertices in the current edge.
          
        u = find(u, parents);//find the parent of u
        v = find(v, parents);//find the parent of v
        if (size.get(u) > size.get(v)) {//merge the larger disjoint set with smaller disjoint set.                                                       
            parents.set(v, u);
            size.set(u, size.get(u) + size.get(v));
        } else {
            parents.set(u, v);
            size.set(v, size.get(v) + size.get(u));
        }

    }
         //////////////////////////////////////////isCyclic method //////////////////////////////////////////////////////////////
     public static boolean isCyclic(int u, int v, LinkedList<Integer> parents) {
        /*
          if the parents of both  vertices of the 
          edge are the same , then this means they are connected 
          to a common vertex. 
           And  if we add this
           edge in the minmum spanning tree then it will create a cycle.
         */
        return find(u, parents) == find(v, parents);
    }
  
        /////////////////////////////////////////////////////Make Graph /////////////////////////////////////////////////    
        public void make_graph(Graph graph) {
            // object of Random class
            Random random = new Random();
           
            for (int j = 0; j < vertices-1; j++) {
                    int weight = random.nextInt(10) + 1;
                    addEdge(j,j+1,weight);
                
            }
            
            // Generate random graph with the remain edges
            int remain = edges- (vertices-1);
            
            for (int j = 0; j < remain; j++) {
                int src = random.nextInt(graph.vertices);
                int destination = random.nextInt(graph.vertices);
                if (destination == src || isConnect(src, destination, graph.adjacencylist)) { // to prevent the duplicate edges  and self loops 
                    j--;
                    continue;
                }
                // Generate random weights in range 0 to 10
                int w = random.nextInt(10) + 1;
                // add  the edge.
                addEdge(src, destination, w);
            }

        }     
        ///////////////////////////////////////////////// is conneted method ////////////////////////////////////////////////////////

        // check if an edge is already existed
        public boolean isConnect(int source, int destination, LinkedList<Edge>[] Edges) {
        for (LinkedList<Edge> j : Edges) {
            for (Edge E : j) {
                if ((E.src == source && E.dest == destination) || (E.src == destination && E.dest == source)) {
                    return true;
                }
            }
        }
        return false;
    }
    }
   
    //////////////////////////////////////////////////MinHeap class ///////////////////////////////////////////////////

    static class Min_heap {
        int size;
        int c_Size;
        HeapN[] min_H;
        int[] ind; 

        public Min_heap(int size) {
            this.size = size;
            min_H = new HeapN[size + 1];
            ind = new int[size];
            min_H[0] = new HeapN();
            min_H[0].k = Integer.MIN_VALUE;
            min_H[0].v = -1;
            c_Size = 0;
        }

        

        public void insert(HeapN N) {
            c_Size++;
            int id = c_Size;
            min_H[id] = N;
            ind[N.v] = id;
            bubbleU(id);
        }

        public void bubbleU(int position) {
            int p_Id = position / 2;
            int c_Id = position;
            while (c_Id > 0 && min_H[p_Id].k > min_H[c_Id].k) {
                HeapN c_Node = min_H[c_Id];
                HeapN p_Node = min_H[p_Id];

                //swap.
                ind[c_Node.v] = p_Id;
                ind[p_Node.v] = c_Id;
                swap(c_Id, p_Id);
                c_Id = p_Id;
                p_Id = p_Id / 2;
            }
        }

        public HeapN eMin() {
            HeapN min = min_H[1];
            HeapN lastNode = min_H[c_Size];
//update the ind[].
            ind[lastNode.v] = 1;
            min_H[1] = lastNode;
            min_H[c_Size] = null;
            s_Down(1);
            c_Size--;//decrese by 1.
            return min;
        }

        public void s_Down(int L) {
            int small = L;
            int left_child_id = 2 * L;
            int right_child_id = 2 * L + 1;
            if (left_child_id < heap_size() && min_H[small].k > min_H[left_child_id].k) {
                small = left_child_id;
            }
            if (right_child_id < heap_size() && min_H[small].k > min_H[right_child_id].k) {
                small = right_child_id;
            }
            if (small != L) {

                HeapN small_N = min_H[small];
                HeapN k_Node = min_H[L];

                //swap.
                ind[small_N.v] = L;
                ind[k_Node.v] = small;
                swap(L, small);
                s_Down(small);
            }
        }

        public void swap(int i, int j) {
            HeapN temp = min_H[i];
            min_H[i] = min_H[j];
            min_H[j] = temp;
        }

        public boolean isEmpty() {
            return c_Size == 0;
        }

        public int heap_size() {
            return c_Size;
        }
    }

    
    
    /////////////////////////////////////////////////////////////main class ///////////////////////////////////////////////////////////////
    
    public static void main(String[] args) {
        int v=0,e=0;
        Scanner in = new Scanner(System.in);
        System.out.println("calculating Runtime of Different Minimum Spanning Tree Algorithms \n\t1- Kruskal's Algorithm & PQ based - Prim's Algorithm "
                + "\n\t2- Min heap based - Prim's Algorithm & PQ based - Prim's Algorithm ");
        System.out.print("> Enter your choice (1 or 2): ");
        int choice = in.nextInt();
        System.out.println("> The cases (where n represents # of vertices and m represents # of edges: ");
        System.out.println(" 1-  n=1,000 - m=10,000");
        System.out.println(" 2-  n=1,000 - m=15,000");
        System.out.println(" 3-  n=1,000 - m=25,000");
        System.out.println(" 4-  n=5,000 - m=15,000");
        System.out.println(" 5-  n=5,000 - m=25,000");
        System.out.println(" 6- n=10,000 - m=15,000");
        System.out.println(" 7  n=10,000 - m=25,000");
        System.out.println(" 8- n=20,000 - m=200,000");
        System.out.println(" 9- n=20,000 - m=300,000");
        System.out.println("10- n=50,000 - m=10,000,000");
        System.out.print("> Enter a case to test: ");
        int input = in.nextInt();
        switch (input) {
            case 1: {
                v=1000; e=10000;
            }
            break;
            case 2: {
                v=1000; e=15000;
            }
            break;
            case 3: {
                v=1000; e=25000;
            }
            break;
            case 4: {
                v=5000; e=15000;
            }
            break;
            case 5: {
                v=5000; e=25000;
            }
            break;
            case 6: {
                v=10000; e=15000;
            }
            break;
            case 7: {
                v=10000; e=25000;
            }
            break;
            case 8: {
                v=20000; e=200000;
            }
            break;
            case 9: {
                v=20000; e=300000;
            }
            break;
            case 10: {
                v=50000; e=1000000;
            }
            break;
            default:
                System.out.println("Wrong choise");
        }
        Graph graph = new Graph(v, e);
        graph.make_graph(graph);
     
        
        switch(choice){
                
            case 1:
                // start time
                double initial_time4 = System.currentTimeMillis();

                graph.kruskalMST(v);//call KruskalMST Method

                //finish time of the algorithm
                double final_time4 = System.currentTimeMillis();

                //print the total time consumed by the algorithm
                System.out.println("Total runtime of Kruskal's Algorithm: " + (final_time4 - initial_time4) + " ms.");

                //////////////////////////////////// the second algorithm /////////////////////////////////
                //start time
                double initial_time3 = System.currentTimeMillis();

                graph.primpq();//call primpq method

                //finish time of the algorithm
                double final_time3 = System.currentTimeMillis();

                //print the total time consumed by the algorithm
                System.out.println("\"Total runtime of Prim's Algorithm (Usin PQ): " + (final_time3 - initial_time3) + " ms.");
                break;

            case 2:

                //start time
                double initial_time = System.currentTimeMillis();

                graph.MinHeap_Based_Prime();//call MinHeap_Based_Prime() method

                //finish time of the algorithm
                double final_time = System.currentTimeMillis();

                //print the total time consumed by the algorithm
                System.out.println("Total runtime of Prim's Algorithm (Usin Min Heap): " + (final_time - initial_time) + " ms.");

                //////////////////////////////////// the second algorithm /////////////////////////////////
                //start time
                double initial_time2 = System.currentTimeMillis();

                graph.primpq();//call primpq()method.

                //finish time of the algorithm
                double final_time2 = System.currentTimeMillis();

                //print the total time consumed by the algorithm
                System.out.println("\"Total runtime of Prim's Algorithm (Usin PQ): " + (final_time2 - initial_time2) + " ms.");

                break;

            default:
                System.out.println("Wrong choise");

        }
                
    }
    
}
    



