import PriorityQueue from './PriorityQueue';

class Graph {
  /**
   * Create a graph.
   * @param {boolean} dir - Specify whether graph is directed or not.
   * @param {boolean} weighted - Specify whether graph is weighted or not.
   * NOTE: Params should be in the form of an object!!
   *       Enter an empty object if configuration is unnecessary.
   */
  constructor({ dir = false, weighted = false }) {
    this.adjacencyList = {};
    this.dir = dir;
    this.weighted = weighted;
  }

  /**
   *
   * @param {number|string} vertex - Vertex to add to graph.
   */
  addVertex(vertex) {
    // Check if vertex is already in graph
    if (!this.adjacencyList[vertex]) this.adjacencyList[vertex] = [];
  }

  /**
   *
   * @param {number|string} src - The source vertex of the edge.
   * @param {number|string} target - The target vertex of the edge.
   */
  addEdge(src, target, weight = null) {
    // Check if the vertices exist in graph
    if (this.adjacencyList[src] && this.adjacencyList[target]) {
      const t = this.weighted
        ? { target, weight }
        : { target, weight: undefined }; // If weighted -> add object with weight
      this.adjacencyList[src].push(t);
      // Add opposite edge if graph is directed
      if (!this.dir) {
        const s = this.weighted
          ? { target: src, weight }
          : { target: src, weight: undefined }; // If weighted -> add object with weight
        this.adjacencyList[target].push(s);
      }
    }
  }

  /**
   *
   * @param {number|string} src - The source vertex of the edge to be removed.
   * @param {number|string} target - The target vertex of the edge to be removed.
   */
  removeEdge(src, target) {
    // Check if vertices exist in graph
    if (this.adjacencyList[src] && this.adjacencyList[target]) {
      this.adjacencyList[src] = this.adjacencyList[src].filter(
        v => v.target !== target
      );
      // Remove opposite edge if graph is directed
      if (!this.dir) {
        this.adjacencyList[target] = this.adjacencyList[target].filter(
          v => v.target !== src
        );
      }
    }
  }

  /**
   *
   * @param {number|string} vertex - The vertex to be removed.
   */
  removeVertex(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex]) {
      // Loop over and remove all edges then delete vertex from graph
      while (this.adjacencyList[vertex].length) {
        const v = this.adjacencyList[vertex].pop();
        this.removeEdge(vertex, v.target);
      }
      delete this.adjacencyList[vertex];
    }
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to get successors.
   * @return {array} Array of successors (undefined if vertex doesn't exist).
   */
  getSucc(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex])
      return this.adjacencyList[vertex].reduce((acc, curr) => {
        acc.push(curr.target);
        return acc;
      }, []);

    // Return undefined if vertex doesn't exist
    return undefined;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to get predecessors.
   * @return {array} Array of predecessors (undefined if vertex doesn't exist).
   */
  getPred(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex]) {
      if (!this.dir)
        return this.adjacencyList[vertex].reduce((acc, curr) => {
          acc.push(curr.target);
          return acc;
        }, []);

      return Object.keys(this.adjacencyList).reduce((acc, curr) => {
        if (this.adjacencyList[curr].map(v => v.target).indexOf(vertex) !== -1)
          acc.push(curr);
        return acc;
      }, []);
    }

    // Return undefined if vertex doesn't exist
    return undefined;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to get neighbors.
   * @return {array} Array of neighbors (undefined if vertex doesn't exist).
   */
  getNeighbors(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex]) {
      if (!this.dir)
        return this.adjacencyList[vertex].reduce((acc, curr) => {
          acc.push(curr.target);
          return acc;
        }, []);
      return [...this.getSucc(vertex), ...this.getPred(vertex)];
    }

    // Return undefined if vertex doesn't exist
    return undefined;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to calculate outdegree.
   * @return {number} Value of outdegree (-1 if vertex doesn't exist).
   */
  outDegree(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex]) return this.adjacencyList[vertex].length;

    // Return -1 if vertex doesn't exist
    return -1;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to calculate indegree.
   * @return {number} Value of indegree (-1 if vertex doesn't exist).
   */
  inDegree(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex]) {
      return Object.keys(this.adjacencyList).reduce((acc, curr) => {
        this.adjacencyList[curr].forEach(
          v => (acc = v.target === vertex ? acc + 1 : acc)
        );
        return acc;
      }, 0);
    }

    // Return -1 if vertex doesn't exist
    return -1;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to calculate degree.
   * @return {number} Value of degree (-1 if vertex doesn't exist).
   */
  degree(vertex) {
    // Check if vertex is in graph
    if (this.adjacencyList[vertex]) {
      if (!this.dir) return this.outDegree(vertex); // degree = outDegree if graph is undirected

      return this.inDegree(vertex) + this.outDegree(vertex);
    }

    // Return -1 if vertex doesn't exist
    return -1;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to check if pendant.
   * @return {boolean} Whether or not the vertex is pendant.
   */
  isPendant(vertex) {
    if (this.degree(vertex) === 1) return true;
    return false;
  }

  /**
   *
   * @param {number|string} vertex - Vertex for which to check if isolated.
   * @return {boolean} Whether or not the vertex is isolated.
   */
  isIsolated(vertex) {
    if (this.degree(vertex) === 0) return true;
    return false;
  }

  /**
   *
   * @param {number|string} start - Vertex to start depth first search.
   * @return {array} Vertices visited in order.
   */
  dfs(start) {
    const result = [];
    const visited = {};
    const adjacencyList = this.adjacencyList;
    let finishOrder = [];

    (function helper(vertex) {
      if (!vertex) return null;
      visited[vertex] = true;
      result.push(vertex);
      adjacencyList[vertex].forEach(v => {
        if (!visited[v.target]) {
          return helper(v.target);
        }
      });
      finishOrder.push(vertex);
    })(start);

    return { result, finishOrder };
  }

  /**
   * @return {boolean} Whether or not the graph is connected.
   */
  isConnected() {
    // Check whether or not you can visit all vertices from first vertex.
    if (
      this.dfs(Object.keys(this.adjacencyList)[0]).result.length ===
      Object.keys(this.adjacencyList).length
    )
      return true;
    return false;
  }

  /**
   * @return {boolean} Whether or not the graph is connected.
   */
  isSimple() {
    let loop = false;
    let multiEdge = false;
    Object.keys(this.adjacencyList).forEach(vertex => {
      const visited = {}; // Object to store "visited" vertices for nested loop.
      this.adjacencyList[vertex].forEach(v => {
        if (visited[v.target]) {
          // Check if already "visited" -> has multiple edges.
          multiEdge = true;
        } else {
          visited[v.target] = 1;
        }
        if (vertex === v.target) loop = true; // If vertex is the same as target -> has loop.
      });
    });

    return !(loop || multiEdge);
  }

  /**
   * @return {boolean} Whether or not the graph is complete.
   */
  isComplete() {
    const graphSize = Object.keys(this.adjacencyList).length;
    const complete = Object.keys(this.adjacencyList).some(vertex => {
      const visited = {}; // Object to store "visited" vertices for nested loop.
      this.adjacencyList[vertex].forEach(v => {
        // Keep track of which vertices have been "visited".
        if (!visited[v.target]) visited[v.target] = 1;
      });
      const visitedLength = Object.keys(visited).length;

      // Array.some() function returns true if condition is meet at least once
      // Therefore if atleast one vertex has less than the graphSize - 1,
      // The graph must be uncomplete
      return !(visitedLength !== graphSize - 1);
    });

    return complete;
  }

  /**
   * @return {boolean} - Whether or not the graph is Eulerien
   */
  isEuler() {
    if (!this.dir) {
      return !Object.keys(this.adjacencyList).some(v => {
        return this.degree(v) % 2 === 1;
      });
    } else {
      return !Object.keys(this.adjacencyList).some(v => {
        return this.inDegree(v) !== this.outDegree(v);
      });
    }
  }

  /**
   * @return {object} - Reversed adjacency list.
   */
  reverseGraph() {
    let reverseAdj = {};
    Object.keys(this.adjacencyList).forEach(vertex => {
      reverseAdj[vertex] = [];
    });

    Object.keys(this.adjacencyList).forEach(vertex => {
      this.adjacencyList[vertex].forEach(v => {
        reverseAdj[v.target] = [
          ...reverseAdj[v.target],
          { target: vertex, weight: v.weight }
        ];
      });
    });
    this.adjacencyList = { ...reverseAdj };

    return reverseAdj;
  }

  /**
   * @return {[array]} - An array of arrays containing all strongly connected components
   */
  stronglyConnected() {
    let finishedStack = [];
    let adjList = { ...this.adjacencyList };

    // Keep track of finish time for each vertex in dfs
    while (Object.keys(adjList).length > 0) {
      const start = Object.keys(adjList)[0];
      const finishOrder = this.dfs(start).finishOrder;
      finishOrder.forEach(v => {
        if (finishedStack.indexOf(v) === -1) finishedStack.push(v);
        delete adjList[v];
      });
    }

    // Reverse the graph and perform dfs again to find SCCs
    this.reverseGraph();
    let visited = {};
    let sccs = [];
    while (finishedStack.length > 0) {
      const start = finishedStack.pop();
      if (!visited[start]) {
        const component = this.dfs(start).finishOrder.filter(v => !visited[v]);
        sccs.push(component);
        component.forEach(v => (visited[v] = true));
      }
    }

    // Reverse the graph back to original form
    this.reverseGraph();
    return sccs;
  }

  /**
   *
   * @param {number|string} start - Beginning vertex
   * @param {number|string} finish - Ending vertex
   * @return {array} - Array containing vertices in shortest path in order
   */
  dijkstra(start, finish) {
    // If graph isn't weighted, can't use dijkstra
    if (!this.weighted) return undefined;

    const nodes = new PriorityQueue();
    const distances = {};
    const previous = {};
    let smallest;
    let path = []; // To return at end

    // Build up initial state
    for (let vertex in this.adjacencyList) {
      if (vertex === start) {
        distances[vertex] = 0;
        nodes.enqueue(vertex, 0);
      } else {
        distances[vertex] = Infinity;
        nodes.enqueue(vertex, Infinity);
      }
      previous[vertex] = null;
    }

    // As long as there is something to visit
    while (nodes.values.length) {
      smallest = nodes.dequeue().val;
      if (smallest === finish) {
        // WE ARE DONE
        // BUILD UP PATH TO RETURN AT END
        while (previous[smallest]) {
          path.push(smallest);
          smallest = previous[smallest];
        }
        break;
      }
      if (smallest || distances[smallest] !== Infinity) {
        for (let neighbor in this.adjacencyList[smallest]) {
          // Find neighboring node
          let nextNode = this.adjacencyList[smallest][neighbor];

          // Calculate new distance to neighboring node
          let candidate = distances[smallest] + nextNode.weight;
          let nextNeighbor = nextNode.target;
          if (candidate < distances[nextNeighbor]) {
            // Updating new smallest distance
            distances[nextNeighbor] = candidate;
            // Updating previous - How we got to neighbor
            previous[nextNeighbor] = smallest;
            // Enqueue in priority queue with new priority
            nodes.enqueue(nextNeighbor, candidate);
          }
        }
      }
    }
    return path.concat(start).reverse();
  }

  /**
   * @param {number|string} startVertex - Vertex to start algorithm from
   * @return {{distances, previousVertices}} - Object containing list of distances and previous vertices
   * (returns undefined if graph isn't weighted)
   */
  bellmanFord(startVertex) {
    // If graph isn't weighted, can't use Bellman-Ford
    if (!this.weighted) return undefined;

    const distances = {};
    const previousVertices = {};

    // Init all distances with infinity assuming that currently we can't reach
    // any of the vertices except start one.
    distances[startVertex] = 0;
    Object.keys(this.adjacencyList).forEach(vertex => {
      previousVertices[vertex] = null;
      if (vertex !== startVertex) {
        distances[vertex] = Infinity;
      }
    });

    // We need (|V| - 1) iterations.
    for (let i = 0; i < Object.keys(this.adjacencyList).length - 1; i++) {
      // During each iteration go through all vertices.
      Object.keys(distances).forEach(vertexKey => {
        // const vertex = this.adjacencyList[vertexKey];

        // Go through all vertex edges.
        this.adjacencyList[vertexKey].forEach(neighbor => {
          // Find out if the distance to the neighbor is shorter in this iteration
          // then in previous one.
          const distanceToVertex = distances[vertexKey];
          const distanceToNeighbor = distanceToVertex + neighbor.weight;
          if (distanceToNeighbor < distances[neighbor.target]) {
            distances[neighbor.target] = distanceToNeighbor;
            previousVertices[neighbor.target] = vertexKey;
          }
        });
      });
    }

    return {
      distances,
      previousVertices
    };
  }

  /**
   * @return {array} - Array containing {src, target, val}, the edges of the minimum spanning tree
   * (returns undefined if graph isn't weighted or isn't connected)
   */
  prim() {
    // If graph isn't weighted or connected, can't use Prim
    if (!this.weighted || !this.isConnected()) return undefined;

    const visited = {};
    const result = [];
    // Add first vertex to visited list
    visited[Object.keys(this.adjacencyList)[0]] = true;
    Object.keys(this.adjacencyList).forEach(vertex => {
      let minWeight = { src: undefined, target: undefined, val: Infinity };
      Object.keys(visited).forEach(v => {
        this.adjacencyList[v].forEach(neighbor => {
          if (neighbor.weight < minWeight.val && !visited[neighbor.target])
            minWeight = {
              src: v,
              target: neighbor.target,
              val: neighbor.weight
            };
        });
      });
      visited[minWeight.target] = true;
      if (minWeight.src) result.push(minWeight);
    });
    return result;
  }
}

export default Graph;
