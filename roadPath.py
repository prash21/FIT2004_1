# PRASHANT MURALI - 29625564
# FIT2004 Assignment 4
import math

# Edge class
class Edge:
    def __init__(self, u, v, w):
        """
        This function is to initialize the edges. It holds the vertex u, adjacent vertex v, its time
        between the two vertices (weight), and a boolean flag to hold if this edge has a toll or not (True of False).
        Time complexity: O(1)
        Space complexity: O(1)
        Error handle: None
        Return: None
        Parameter: vertex u, vertex v, and time w(represented as 'weight')
        Precondition: u represented as an integer, v represented as an integer, time w in float
        """
        self.u = u
        self.v = v
        self.weight=w
        self.toll=False


# Vertex class
class Vertex:
    """
    This function is to initialize the vertices. It holds the vertex u, and a boolean flag to hold if
    this vertex has a camera or not (True or False). It also has another boolean flag to hold if this
    vertex has a service or not (True or False).
    Time complexity: O(1)
    Space complexity: O(1)
    Error handle: None
    Return: None
    Parameter: vertex u
    Precondition: u represented as an integer
    """
    def __init__(self, u):
        self.u=u
        self.camera=False
        self.service=False


# Graph class
class Graph:
    def __init__(self):
        """
        This function is to initialize the graph. It would be used to represent a road network.
        Time complexity: O(1)
        Space complexity: O(E+V)
        Error handle: None
        Return: None
        Parameter: None
        Precondition: None
        """
        self.adjacency_list=[]
        self.N=0
        self.vertex_list=[]
        self.service_data=None
        self.data_list=None

    def buildGraph(self, filename_roads):
        """
        This function populates the Graph object with vertices and edges according to the edge list present
        in the file provided.
        Time complexity: O(E+V)
        Space complexity: O(E+V)
        Error handle: None (filename errors are handled in the main block)
        Return: None
        Parameter: filename_roads which is the name of the file which contains the vertex u, v and time w
        Precondition: vertex ID in file is continious, starting from vertex 0 up to vertex k where k is the
                      largest vertex in the edge list sorted within the input file.
        """
        file = open(filename_roads, "r")

        # Read the file.
        file = file.read()

        # Put all the data from the input file into data_list.
        self.data_list=self.splitNewLine(file)

        # Get max
        i=0
        max=int(self.data_list[0])
        while i<len(self.data_list):
            if int(self.data_list[i])>max:
                max=int(self.data_list[i])
            i+=3

        # Set self.N to max+1 to denote the number of vertices in the graph.
        self.N=max+1

        # Expand the adjacency list and vertex list accordingly.
        for i in range(max+1):
            self.adjacency_list.append([])
            self.vertex_list.append([])

        # Place the data into their respective classes and lists.
        index=0
        while index<len(self.data_list):
            i=0
            u=(int(self.data_list[index]))
            v=(int(self.data_list[index+1]))
            w=(float(self.data_list[index+2]))
            self.adjacency_list[u].append(Edge(u,v,w))
            index+=3

        # Add the Vertex class data to the vertex_list
        index=0
        counter=0
        vertex_list2 = []
        while index<len(self.data_list):
            if counter<2:
                vertex_list2.append(self.data_list[index])
                counter+=1
            else:
                counter=0
            index+=1

        for item in vertex_list2:
            self.vertex_list[int(item)]=(Vertex(int(item)))

    def buildReverseGraph(self, my_data_list):
        """
        Note: This function is used for task 3.
        This function populates the Graph object with vertices and edges according to the edge list present
        data provided. The vertex v and vertex u are swapped before adding int the graph, to
        make this the reverse graph.
        Time complexity: O(E+V)
        Space complexity: O(E+V)
        Error handle: None (filename errors are handled in the main block)
        Return: None
        Parameter: my_data_list is the list of data of vertices that were put in to the file earlier.
        Precondition: vertex ID in file is continuous, starting from vertex 0 up to vertex k where k is the
                      largest vertex in the edge list sorted within the input file.
        """
        data_list=my_data_list

        # Get max
        i = 0
        max = int(data_list[0])
        while i < len(data_list):
            if int(data_list[i]) > max:
                max = int(data_list[i])
            i += 3

        # Set self.N to max+1 to denote the number of vertices in the graph.
        self.N = max + 1

        # Expand the adjacency list and vertex list accordingly.
        for i in range(max + 1):
            self.adjacency_list.append([])
            self.vertex_list.append([])

        # Place the data into their respective classes and lists.
        index = 0
        while index < len(data_list):
            i = 0
            v = (int(data_list[index]))
            u = (int(data_list[index + 1]))
            w = (float(data_list[index + 2]))
            self.adjacency_list[u].append(Edge(u, v, w))
            index += 3

        # Add the Vertex class data to the vertex_list
        index = 0
        counter = 0
        vertex_list2 = []
        while index < len(data_list):
            if counter < 2:
                vertex_list2.append(data_list[index])
                counter += 1
            else:
                counter = 0
            index += 1

        for item in vertex_list2:
            self.vertex_list[int(item)] = (Vertex(int(item)))


    def splitNewLine(self,file):
        """
        This function splits items in the file (between newlines).
        Time complexity: O(E+V)
        Space complexity: O(E+V)
        Error handle: None
        Return: A list with each item in it.
        Parameter: file, which is the input file given.
        Precondition: The parameter must be a file with items that need to be split between newlines.
        """
        data_list=[]
        temp=""
        for index in range(len(file)):
            if file[index] != " ":
                if file[index] != "\n":
                    temp+=file[index]
                    if index == (len(file) - 1):
                        data_list.append(temp)
                else:
                    data_list.append(temp)
                    temp = ""
            else:
                data_list.append(temp)
                temp=""

        return data_list

    # TASK 1
    def quickestPath(self, source, target):
        """
        This function finds the quickest path from a given starting location u, to destination v.
        Time complexity: O(ElogV)
        Space complexity: O(E+V)
        Error handle: None
        Return: List containing all nodes in order of quickest path traversal from source to target, and
                the time taken, returned as a tuple.
        Parameter: source which denotes the starting point of the travel, and target which denotes the destination
                   point of travel.
        Precondition: Both source and target are valid vertices in the graph.
        """
        # Initialize the min heap
        min_heap = minHeap()

        # Get the number of vertices in graph
        num_v=self.N
        # Used to get time (weight of the edge)
        time = []
        # Used to hold the path
        path=[]

        # Every node of min heap should contain vertex number and distance value of the vertex
        for v in range(num_v):
            path.append(None)
            time.append(math.inf)
            min_heap.array.append([v, time[v]])
            min_heap.pos.append(v)

        # Make the source vertex as the root. Time value for root would be 0.
        min_heap.pos[source] = source
        time[source] = 0
        min_heap.decreaseKey(source, time[source])
        min_heap.size = num_v

        # Now, extract the vertex with minimum time value node from the heap and for every
        # adjacent vertex v of u, check if v is in the heap.
        # If v is in heap and its time value is more than the time taken (weight) of u-v plus time
        # value of u, then the time value of v is updated.
        while min_heap.size != 0:

            # Extracting vertex with minimum time value
            newHeapNode=min_heap.extractMin()
            u=newHeapNode[0]

            # Traverses through all adjacent vertices of
            # u and updates their time values accordingly.
            for edge in self.adjacency_list[u]:
                v = edge.v

                # Given that vertex u is now the target, the loop can be broken.
                if u == target:
                    break

                # Updating time
                # Time value in min heap is also updated.
                if min_heap.pos[v]<min_heap.size:
                    if time[u] != math.inf:
                        if edge.weight + time[u] < time[v]:
                            time[v] = edge.weight + time[u]
                            # Keep track of path
                            path[v]=u
                            min_heap.decreaseKey(v, time[v])

        # time list now has all shortest time (weight) from source.
        # Extract the shortest time (weight) from source to target.
        shortest_time=(time[target])

        # If the shortest time is infinity, it simply means a path does not exist
        # between source and target.
        if shortest_time==math.inf:
            return ([[],-1])

        # Now extract the path from the path list.
        quickest_path = self.getPath(path, target)

        # Returned as a tuple.
        return quickest_path,shortest_time

    # TASK 2
    def augmentGraph(self, filename_camera, filename_toll):
        """
        This function is used to add in the additional elements into the Graph data structure for Task2.
        Time complexity: O(E+V)
        Space complexity: O(E+V)
        Error handle: None
        Return: None
        Parameter: filename_camera which denotes filename of the text file containing the list of red light cameras;
                   and filename_toll which denotes the filename of the text file containing the list of tolls.
        Precondition: The vertices listed in the file are all valid vertices.
        """

        # Open and read the files.
        camera_file = open(filename_camera, "r")
        toll_file = open(filename_toll, "r")

        camera_file = camera_file.read()
        toll_file = toll_file.read()

        # Remove new lines from the file and place it into a list.
        camera_data=self.splitNewLine(camera_file)

        # Remove new lines from the file and place it into a list.
        toll_data = self.splitNewLine(toll_file)

        # Alter the 'self.camera' value from the vertex class in the vertex_list,
        # so that the vertices that were given in the camera file will have the
        # 'self.camera' variable set to True, indicating that is has a red light camera.
        for item in camera_data:
            self.vertex_list[int(item)].camera=True

        # Similarly done as above.
        # For all the edges in the adjacency_list, mark the edge's 'self.toll' to True
        # for the given u,v pairs in the toll file.
        i=0
        while i<len(toll_data):
            for item in self.adjacency_list[int(toll_data[i])]:
                if item.v==int(toll_data[i+1]):
                    item.toll=True
            i+=2

    def quickestSafePath(self, source, target):
        """
        Finds the quickest safe path such that none of the vertices within the path contains
        a red light camera anf none of the edges are toll roads.
        Time complexity: O(ElogV)
        Space complexity: O(E+V)
        Error handle: None
        Return: List containing all nodes in order of quickest safe path traversal from source to target, and
                the time taken, returned as a tuple.
        Parameter: source which denotes the starting point of the travel, and target which denotes the destination
                   point of travel.
        Precondition: Both source and target are valid vertices in the graph.
        """
        # Initialize the min heap
        min_heap = minHeap()

        # Get the number of vertices in graph
        num_v = self.N
        # Used to get time (weight of the edge)
        time = []
        # Used to hold the path
        path = []

        # Every node of min heap should contain vertex number and distance value of the vertex
        for v in range(num_v):
            path.append(None)
            time.append(math.inf)
            min_heap.array.append([v, time[v]])
            min_heap.pos.append(v)

        # Make the source vertex as the root. Time value for root would be 0.
        min_heap.pos[source] = source
        time[source] = 0
        min_heap.decreaseKey(source, time[source])
        min_heap.size = num_v

        # Now, extract the vertex with minimum time value node from the heap and for every
        # adjacent vertex v of u, check if v is in the heap.
        # If v is in heap and its time value is more than the time taken (weight) of u-v plus time
        # value of u, and if the vertex does not have a red light camera and its edge is not a toll,
        # then the time value of v is updated.
        while min_heap.size != 0:

            # Extract the vertex with minimum time value
            newHeapNode=min_heap.extractMin()
            u=newHeapNode[0]

            # Traverses through all adjacent vertices of
            # u and updates their time values accordingly.
            for edge in self.adjacency_list[u]:
                v = edge.v

                # Given that vertex u is now the target, the loop can be broken.
                if u == target:
                    break

                # Updating time
                # Time value in min heap is also updated.
                if min_heap.pos[v]<min_heap.size:
                    if time[u] != math.inf:
                        if edge.weight + time[u] < time[v]:
                            # Make sure that the vertex does not have a red light camera.
                            # Also make sure the edge does not have a toll.
                            if edge.toll != True and (self.vertex_list[u]).camera != True:
                                time[v] = edge.weight + time[u]
                                # Keep track of path
                                path[v] = u
                                min_heap.decreaseKey(v, time[v])

        # time list now has all shortest time (weight) from source.
        # Extract the shortest time (weight) from source to target.
        shortest_time = (time[target])

        # If the shortest time is infinity, it simply means a path does not exist
        # between source and target.
        if shortest_time == math.inf:
            return ([[], -1])


        # Now extract the path from the path list.
        quickest_path = self.getPath(path,target)

        # Returned as a tuple.
        return quickest_path, shortest_time

    # TASK 3
    def addService(self, filename_service):
        """
        This function is used to add in the additional service elements into the Graph data structure for Task3.
        Time complexity: O(V)
        Space complexity: O(V)
        Error handle: None
        Return: None
        Parameter: filename_service which denotes filename of the text file containing the list of vertices that
                   are service points.
        Precondition: The vertices listed in the file are all valid vertices.
        """

        # Open and read the files.
        service_file = open(filename_service, "r")

        service_file = service_file.read()

        # Remove new lines from the file and place it into a list.
        self.service_data = self.splitNewLine(service_file)

        # Alter the 'self.service' value from the vertex class in the vertex_list,
        # so that the vertices that were given in the service file will have the
        # 'self.service' variable set to True, indicating that it is a service point.
        for item in self.service_data:
            self.vertex_list[int(item)].service = True


    def getPath(self, path, target):
        """
        This function extracts the path from the parent path list produced from the Dijkstra's algorithm.
        Time complexity: O(V)
        Space complexity: O(V)
        Error handle: None
        Return: None
        Parameter: path which is the list of paths produced by the Dijkstra's algorithm, and target
                   which is the target node.
        Precondition: None
        """
        quickest_path = []
        quickest_path.append(target)
        i = target

        # Get the path
        bflag = True
        while bflag == True:
            if path[i] != target and path[i] != None:
                quickest_path.append(path[i])
                i = path[i]
            else:
                bflag = False

        # Reverse the path to get final path - as path was read from
        # backwards just now.
        quickest_path = (quickest_path[::-1])
        return quickest_path

    def quickestDetourPath(self, source, target):
        """
        This function finds the path from the source vertex to the target vertex which passes through
        at least one of the service vertices.
        Time complexity: O(ELogV)
        Space complexity: O(E+V)
        Error handle: None
        Return: None
        Parameter: List containing all nodes in order of quickest safe path traversal from source to target that
                   passes through at least one service point, and the time taken, returned as a tuple.
        Precondition: Both source and target are valid vertices in the graph.
        """
        # Run dijkstra's from source to target.
        path1,time1,num1=(self.task_c_quickest_path(source, target))

        # Get the shortest path from source to target.
        path=(self.getPath(path1, target))

        # If the shortest path already contains a service vertex, return the path
        # and shortest time taken.
        myflag=False
        for vertex in path:
            if self.vertex_list[vertex].service==True:
                myflag=True

        # return the path and shortest time taken if service point was found in shortest path.
        if myflag==True:
            return path,time1[target]

        # Else look for the detour path.
        else:
            # Create another graph to store the vertices but its direction is flipped
            # this time.
            reverseGraph=Graph()
            reverseGraph.buildReverseGraph(self.data_list)

            # Run dijkstra's from the target to source now.
            path2,time2,num2=(reverseGraph.task_c_quickest_path(target,source))

            current_total_time=0
            new_target=None

            # Look through all the service points.
            for i in self.service_data:
                i=int(i)
                # Check if the service point is set to True in the vertex_list.
                if self.vertex_list[i].service == True:
                        # If so, check if at service point i, both forward and
                        # reverse graph have valid paths.
                        if time2[i]!=math.inf and time1[i]!=math.inf:
                            # If so, check if the total time take  is less than the
                            # previous path with service point (if any).
                            if time1[i]+time2[i]>current_total_time:
                                # If so, current time and new_target will be updated.
                                current_total_time=time1[i]+time2[i]
                                new_target = i

            # If the new target remains as None, it means that theres no possible
            # detour path that passes through the service points.
            if new_target == None:
                return [[], -1]

            # If there is, get the paths from source to new_target, and target to new_target.
            else:
                path3, time3, num3 = (self.task_c_quickest_path(source, new_target))
                path_3 = (self.getPath(path3, new_target))

                path4, time4, num4 = (reverseGraph.task_c_quickest_path(target, new_target))
                path_4 = (self.getPath(path4, new_target))

                # Join the paths together.
                i = -2
                while i >= (-(len(path_4))):
                    path_3.append(path_4[i])
                    i -= 1

                # Path and total time taken is returned.
                return path_3, current_total_time



    def task_c_quickest_path(self, source, target):
        """
        Note: This is essentially the implementation of Dijkstra from the quickestPath function
               but has its returned values changed a little for task c.
        Finds the quickest safe path such that none of the vertices within the path contains
        a red light camera anf none of the edges are toll roads.
        Time complexity: O(ElogV)
        Space complexity: O(E+V)
        Error handle: None
        Return: Path parent list, time list and number of vertices returned as a tuple.
        Parameter: source which denotes the starting point of the travel, and target which denotes the destination
                   point of travel.
        Precondition: Both source and target are valid vertices in the graph.
        """
        # Initialize the min heap
        min_heap = minHeap()

        # Get the number of vertices in graph
        num_v = self.N
        # Used to get time (weight of the edge)
        time = []
        # Used to hold the path
        path = []

        # Every node of min heap should contain vertex number and distance value of the vertex
        for v in range(num_v):
            path.append(None)
            time.append(math.inf)
            min_heap.array.append([v, time[v]])
            min_heap.pos.append(v)

        # Make the source vertex as the root. Time value for root would be 0.
        min_heap.pos[source] = source
        time[source] = 0
        min_heap.decreaseKey(source, time[source])
        min_heap.size = num_v

        # Now, extract the vertex with minimum time value node from the heap and for every
        # adjacent vertex v of u, check if v is in the heap.
        # If v is in heap and its time value is more than the time taken (weight) of u-v plus time
        # value of u, then the time value of v is updated.
        # newflag=False
        while min_heap.size != 0:

            # Extracting vertex with minimum time value
            newHeapNode = min_heap.extractMin()
            u = newHeapNode[0]

            # Traverses through all adjacent vertices of
            # u and updates their time values accordingly.
            for edge in self.adjacency_list[u]:
                v = edge.v

                # Given that vertex u is now the target, the loop can be broken.
                if u == target:
                    break

                # Updating time
                # Time value in min heap is also updated.
                if min_heap.pos[v] < min_heap.size:
                    if time[u] != math.inf:
                        if edge.weight + time[u] < time[v]:
                            time[v] = edge.weight + time[u]
                            # Keep track of path
                            path[v] = u
                            min_heap.decreaseKey(v, time[v])

        return (path, time, num_v)


# This is an implementation of a min heap and it used when
# implementing dijkstra's algorithm.
class minHeap():

    def __init__(self):
        """
        Initialize the min heap class.
        Time complexity: O(E+V)
        Space complexity: O(E+V)
        Error handle: None
        Return: None
        Parameter: filename_camera which denotes filename of the text file containing the list of red light cameras;
                   and filename_toll which denotes the filename of the text file containing the list of tolls.
        Precondition: The vertices listed in the file are all valid vertices.
        """
        self.array = []
        self.size = 0
        self.pos = []

    def swapMinHeapNode(self, x, y):
        """
        Swaps two nodes within the heap.
        Time complexity: O(1)
        Space complexity: O(1)
        Error handle: None
        Return: None
        Parameter: x and y, which are the two nodes that need to be swapped.
        Precondition: None
        """
        t = self.array[x]
        self.array[x] = self.array[y]
        self.array[y] = t

    def minHeapify(self, index):
        """
        Heapify for index and updates the position of nodes after they're swapped.
        Time complexity: O(logV)
        Space complexity: O(1)
        Error handle: None
        Return: None
        Parameter: Index that needs to be heapified.
        Precondition: None
        """
        # Initialize the right, left and smallest value.
        right = 2 * index + 2
        left = 2 * index + 1
        smallest = index

        if right < self.size and self.array[right][1] \
                < self.array[smallest][1]:
            smallest = right

        if left < self.size and self.array[left][1] \
                < self.array[smallest][1]:
            smallest = left

        # Swap nodes given that index is not the smallest.
        if smallest != index:

            # Swap the positions before swapping the nodes.
            self.pos[self.array[smallest][0]] = index
            self.pos[self.array[index][0]] = smallest
            self.swapMinHeapNode(smallest, index)
            self.minHeapify(smallest)

    def extractMin(self):
        """
        This function extracts minimum node.
        Time complexity: O(logV)
        Space complexity: O(1)
        Error handle: None
        Return: None
        Parameter: None
        Precondition: None
        """
        # Return NULL wif heap is empty
        if self.size==0:
            return

        # Store the root node
        root=self.array[0]

        # Replace root node with last node and update its position.
        lastNode=self.array[self.size-1]
        self.array[0]=lastNode
        self.pos[lastNode[0]]=0
        self.pos[root[0]]=self.size-1

        # Heapify again after reducing the size of the heap.
        self.size-=1
        self.minHeapify(0)

        # Return the root.
        return root

    def decreaseKey(self, v, time):
        """
        Traverses through the heap while it is not heapified, and swaos the node
        with its parent.
        Time complexity: O(LogV)
        Space complexity: O(1)
        Error handle: None
        Return: None
        Parameter: v is the vertex and time is the updated time.
        Precondition: None
        """
        # Get index of v and update its node with the new time value.
        i = self.pos[v]
        self.array[i][1] = time

        # Traverse the heap while not heapified.
        while i > 0 and self.array[i][1] < self.array[(i - 1) // 2][1]:
            # Swap this node with its parent and move to parent index.
            self.pos[self.array[i][0]] = (i - 1) // 2
            self.pos[self.array[(i - 1) // 2][0]] = i
            self.swapMinHeapNode(i, (i - 1) // 2)
            i = (i - 1) // 2


if __name__ == "__main__":

    graph=Graph()

    line="---------------------------------------------------------------------"

    print(line)

    bflag = False
    while bflag == False:
        try:
            graphFileName = input("Enter the file name for the graph : ")
            openTestfile = open(graphFileName, "r")
            bflag = True
        except NameError:
            bflag = False
        except TypeError:
            bflag = False
        except FileNotFoundError:
            bflag = False
        except IOError:
            bflag = False

    bflag = False
    while bflag == False:
        try:
            cameraFileName = input("Enter the file name for camera nodes : ")
            openTestfile = open(cameraFileName, "r")
            bflag = True
        except NameError:
            bflag = False
        except TypeError:
            bflag = False
        except FileNotFoundError:
            bflag = False
        except IOError:
            bflag = False

    bflag = False
    while bflag == False:
        try:
            tollFileName = input("Enter the file name for the toll roads : ")
            openTestfile = open(tollFileName, "r")
            bflag = True
        except NameError:
            bflag = False
        except TypeError:
            bflag = False
        except FileNotFoundError:
            bflag = False
        except IOError:
            bflag = False

    bflag = False
    while bflag == False:
        try:
            serviceFileName = input("Enter the file name for the service nodes : ")
            openTestfile = open(serviceFileName, "r")
            bflag = True
        except NameError:
            bflag = False
        except TypeError:
            bflag = False
        except FileNotFoundError:
            bflag = False
        except IOError:
            bflag = False

    print(line)

    bflag = False
    while bflag == False:
        try:
            source_node = int(input("Source node: "))
            bflag = True
        except TypeError:
            bflag = False
        except ValueError:
            bflag = False
        except IOError:
            bflag = False

    bflag = False
    while bflag == False:
        try:
            target_node = int(input("Sink node: "))
            bflag = True
        except TypeError:
            bflag = False
        except ValueError:
            bflag = False
        except IOError:
            bflag = False

    print(line)

    # For quickest path
    print("Quickest path:")
    graph.buildGraph(graphFileName)
    shortest_path,quickest_time=graph.quickestPath(source_node,target_node)

    if shortest_path!=[]:
        output_str=""
        for vertex in shortest_path:
            output_str+=(str(vertex))
            if vertex!=shortest_path[-1]:
                output_str+=" --> "
        print(output_str)
        print("Time: "+str(quickest_time)+" minute(s)")
    else:
        print("No path exists")
        print("Time: 0 minute(s)")

    print(line)


    # For quickest safe path
    print("Safe quickest path:")
    graph.augmentGraph(cameraFileName,tollFileName)
    shortest_path,quickest_time=graph.quickestSafePath(source_node,target_node)

    if shortest_path!=[]:
        output_str = ""
        for vertex in shortest_path:
            output_str += (str(vertex))
            if vertex != shortest_path[-1]:
                output_str += " --> "
        print(output_str)
        print("Time: " + str(quickest_time) + " minute(s)")
    else:
        print("No path exists")
        print("Time: 0 minute(s)")

    print(line)


    # For quickest detour path
    print("Quickest detour path:")
    graph.addService(serviceFileName)
    shortest_path,quickest_time=((graph.quickestDetourPath(source_node,target_node)))

    if shortest_path!=[]:
        output_str = ""
        for vertex in shortest_path:
            output_str += (str(vertex))
            if vertex != shortest_path[-1]:
                output_str += " --> "
        print(output_str)
        print("Time: " + str(quickest_time) + " minute(s)")
    else:
        print("No path exists")
        print("Time: 0 minute(s)")

    print(line)
    print("Program end")