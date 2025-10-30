// edge_dual_normal_geometry.h
// Created for handling edges with dual normals
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "geometrycentral/utilities/vector3.h"

namespace geometrycentral {

class EdgeDualNormalGeometry {
public:
    // Constructor
    EdgeDualNormalGeometry();
    EdgeDualNormalGeometry(const std::vector<Vector3>& vertices, 
                          const std::vector<std::pair<size_t, size_t>>& edges,
                          const std::vector<Vector3>& normals1,
                          const std::vector<Vector3>& normals2);

    // Destructor
    virtual ~EdgeDualNormalGeometry() = default;

    // Data accessors
    const std::vector<Vector3>& getVertices() const { return V; }
    const std::vector<std::pair<size_t, size_t>>& getEdges() const { return E; }
    const std::vector<Vector3>& getNormals1() const { return N1; }
    const std::vector<Vector3>& getNormals2() const { return N2; }

    // Data setters
    void setVertices(const std::vector<Vector3>& vertices) { V = vertices; }
    void setEdges(const std::vector<std::pair<size_t, size_t>>& edges) { E = edges; }
    void setNormals1(const std::vector<Vector3>& normals1) { N1 = normals1; }
    void setNormals2(const std::vector<Vector3>& normals2) { N2 = normals2; }

    // Utility functions
    size_t nVertices() const { return V.size(); }
    size_t nEdges() const { return E.size(); }
    bool isValid() const;

    // Get normal pair for specific edge
    std::pair<Vector3, Vector3> getEdgeNormals(size_t edgeIdx) const;

    // Print summary
    void printSummary() const;

private:
    std::vector<Vector3> V;  // Vertices
    std::vector<std::pair<size_t, size_t>> E;  // Edges (pairs of vertex indices)
    std::vector<Vector3> N1; // First normal for each edge
    std::vector<Vector3> N2; // Second normal for each edge
};

// Function to read the Python-generated OBJ file
bool readEdgeDualNormal(const std::string& filename, EdgeDualNormalGeometry& geometry);

// Implementation

inline EdgeDualNormalGeometry::EdgeDualNormalGeometry() {}

inline EdgeDualNormalGeometry::EdgeDualNormalGeometry(
    const std::vector<Vector3>& vertices,
    const std::vector<std::pair<size_t, size_t>>& edges,
    const std::vector<Vector3>& normals1,
    const std::vector<Vector3>& normals2)
    : V(vertices), E(edges), N1(normals1), N2(normals2) {}

inline bool EdgeDualNormalGeometry::isValid() const {
    return (E.size() == N1.size()) && (E.size() == N2.size()) && !V.empty();
}

inline std::pair<Vector3, Vector3> EdgeDualNormalGeometry::getEdgeNormals(size_t edgeIdx) const {
    if (edgeIdx >= E.size()) {
        throw std::out_of_range("Edge index out of range");
    }
    return std::make_pair(N1[edgeIdx], N2[edgeIdx]);
}

inline void EdgeDualNormalGeometry::printSummary() const {
    std::cout << "EdgeDualNormalGeometry Summary:" << std::endl;
    std::cout << "- " << nVertices() << " vertices" << std::endl;
    std::cout << "- " << nEdges() << " edges" << std::endl;
    std::cout << "- " << N1.size() + N2.size() << " normal vectors (2 per edge)" << std::endl;
    std::cout << "- Valid: " << (isValid() ? "Yes" : "No") << std::endl;
}

// Reader function implementation
inline bool readEdgeDualNormal(const std::string& filename, EdgeDualNormalGeometry& geometry) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::vector<Vector3> vertices;
    std::vector<std::pair<size_t, size_t>> edges;
    std::vector<Vector3> normals;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string prefix;
        iss >> prefix;

        if (prefix == "v") {
            // Read vertex: v x y z
            double x, y, z;
            if (iss >> x >> y >> z) {
                vertices.push_back(Vector3{x, y, z});
            }
        }
        else if (prefix == "l") {
            // Read edge: l i j (1-based in file, convert to 0-based)
            size_t i, j;
            if (iss >> i >> j) {
                edges.emplace_back(i - 1, j - 1); // Convert to 0-based
            }
        }
        else if (prefix == "vn") {
            // Read normal: vn nx ny nz
            double nx, ny, nz;
            if (iss >> nx >> ny >> nz) {
                normals.push_back(Vector3{nx, ny, nz});
            }
        }
    }

    file.close();

    // Validate data consistency
    if (edges.size() * 2 != normals.size()) {
        std::cerr << "Error: Expected " << edges.size() * 2 
                  << " normals for " << edges.size() << " edges, but got " 
                  << normals.size() << std::endl;
        return false;
    }

    // Split normals into two arrays (first and second normal for each edge)
    std::vector<Vector3> normals1, normals2;
    for (size_t i = 0; i < edges.size(); ++i) {
        normals1.push_back(normals[2 * i]);     // First normal
        normals2.push_back(normals[2 * i + 1]); // Second normal
    }

    // Set the geometry data
    geometry.setVertices(vertices);
    geometry.setEdges(edges);
    geometry.setNormals1(normals1);
    geometry.setNormals2(normals2);

    std::cout << "Successfully read from " << filename << ":" << std::endl;
    geometry.printSummary();

    return true;
}


inline bool resampleEdgeDualNormalGeometry(
    const EdgeDualNormalGeometry& sourceGeometry,
    EdgeDualNormalGeometry& targetGeometry,
    float targetEdgeLength = 1.0f) {
    
    const auto& sourceVertices = sourceGeometry.getVertices();
    const auto& sourceEdges = sourceGeometry.getEdges();
    const auto& sourceNormals1 = sourceGeometry.getNormals1();
    const auto& sourceNormals2 = sourceGeometry.getNormals2();
    
    if (sourceEdges.empty() || sourceVertices.empty()) {
        std::cerr << "Error: Source geometry is empty" << std::endl;
        return false;
    }
    
    if (targetEdgeLength <= 0.0f) {
        std::cerr << "Error: Target edge length must be positive" << std::endl;
        return false;
    }
    
    std::vector<Vector3> newVertices = sourceVertices; // Start with original vertices
    std::vector<std::pair<size_t, size_t>> newEdges;
    std::vector<Vector3> newNormals1;
    std::vector<Vector3> newNormals2;
    
    // Helper function to calculate distance between two vertices
    auto calculateDistance = [](const Vector3& a, const Vector3& b) -> float {
        float dx = b.x - a.x;
        float dy = b.y - a.y;
        float dz = b.z - a.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    };
    
    for (size_t edgeIdx = 0; edgeIdx < sourceEdges.size(); ++edgeIdx) {
        const auto& edge = sourceEdges[edgeIdx];
        size_t startVertexIdx = edge.first;
        size_t endVertexIdx = edge.second;
        
        const Vector3& startVertex = sourceVertices[startVertexIdx];
        const Vector3& endVertex = sourceVertices[endVertexIdx];
        const Vector3& normal1 = sourceNormals1[edgeIdx];
        const Vector3& normal2 = sourceNormals2[edgeIdx];
        
        // Calculate current edge length
        float currentEdgeLength = calculateDistance(startVertex, endVertex);
        
        if (currentEdgeLength < 1e-6f) {
            // Skip degenerate edges (zero length)
            std::cerr << "Warning: Skipping degenerate edge " << edgeIdx << std::endl;
            continue;
        }
        
        // Calculate how many segments we need
        int numSegments = static_cast<int>(std::ceil(currentEdgeLength / targetEdgeLength));
        numSegments = std::max(1, numSegments); // At least 1 segment
        
        // Calculate actual segment length (might be slightly less than target)
        float actualSegmentLength = currentEdgeLength / numSegments;
        
        // Create subdivided segments
        size_t currentVertexIdx = startVertexIdx;
        
        for (int i = 0; i < numSegments; ++i) {
            size_t nextVertexIdx;
            
            if (i == numSegments - 1) {
                // Last segment connects to the original end vertex
                nextVertexIdx = endVertexIdx;
            } else {
                // Create intermediate vertex
                float t = static_cast<float>(i + 1) / numSegments;
                Vector3 interpolatedVertex = {
                    startVertex.x + t * (endVertex.x - startVertex.x),
                    startVertex.y + t * (endVertex.y - startVertex.y),
                    startVertex.z + t * (endVertex.z - startVertex.z)
                };
                
                newVertices.push_back(interpolatedVertex);
                nextVertexIdx = newVertices.size() - 1;
            }
            
            // Add the new edge segment
            newEdges.emplace_back(currentVertexIdx, nextVertexIdx);
            
            // Copy normals to all segments of this edge
            newNormals1.push_back(normal1);
            newNormals2.push_back(normal2);
            
            currentVertexIdx = nextVertexIdx;
        }
    }
    
    // Set the resampled geometry data
    targetGeometry.setVertices(newVertices);
    targetGeometry.setEdges(newEdges);
    targetGeometry.setNormals1(newNormals1);
    targetGeometry.setNormals2(newNormals2);
    
    std::cout << "Resampling complete - Target edge length: " << targetEdgeLength << std::endl;
    std::cout << "Original: " << sourceVertices.size() << " vertices, "
              << sourceEdges.size() << " edges" << std::endl;
    std::cout << "Resampled: " << newVertices.size() << " vertices, "
              << newEdges.size() << " edges" << std::endl;
    
    return true;
}



inline bool snapNormals(EdgeDualNormalGeometry& geometry,
                        float angleTolerance = 15.0f) {
    
    // Get copies of the normal arrays (const methods return const references)
    std::vector<Vector3> normals1 = geometry.getNormals1();
    std::vector<Vector3> normals2 = geometry.getNormals2();
    
    if (normals1.size() != normals2.size()) {
        std::cerr << "Error: Normal arrays must have the same size" << std::endl;
        return false;
    }
    
    if (angleTolerance <= 0.0f || angleTolerance >= 180.0f) {
        std::cerr << "Error: Angle tolerance must be between 0 and 180 degrees" << std::endl;
        return false;
    }
    
    // Convert angle tolerance to cosine for efficient comparison
    float cosAngleTolerance = std::cos(angleTolerance * M_PI / 180.0f);
    
    size_t snappedCount = 0;
    
    for (size_t i = 0; i < normals1.size(); ++i) {
        Vector3 n1 = normals1[i];
        Vector3 n2 = normals2[i];
        
        // n1 and n2 are already normalized unit vectors
        
        // Calculate dot product (cosine of angle between normals)
        float dotProduct = dot(n1, n2);
        
        // Clamp dot product to [-1, 1] to handle numerical errors
        dotProduct = std::max(-1.0f, std::min(1.0f, dotProduct));
        
        // Check if normals are close (angle is small)
        // We use abs(dotProduct) to handle both parallel and anti-parallel cases
        if (std::abs(dotProduct) >= cosAngleTolerance) {
            // Normals are close, compute average
            Vector3 averageNormal = normalize((n1 + n2) / 2.0);
            
            // Set both normals to the average
            normals1[i] = averageNormal;
            normals2[i] = averageNormal;
            
            snappedCount++;
        }
        // If normals are not close enough, keep them as they are (no action needed)
    }
    
    // Set the modified arrays back to the geometry using setter methods
    geometry.setNormals1(normals1);
    geometry.setNormals2(normals2);
    
    std::cout << "Normal snapping complete:" << std::endl;
    std::cout << "- Total edges: " << normals1.size() << std::endl;
    std::cout << "- Snapped normals: " << snappedCount << std::endl;
    std::cout << "- Angle tolerance: " << angleTolerance << " degrees" << std::endl;
    
    return true;
}


inline Vector3 rotateAroundAxis(const Vector3& v, const Vector3& axis, float theta) {
    Vector3 a = normalize(axis);
    float c = cos(theta);
    float s = sin(theta);
    return v * c + cross(a, v) * s + a * dot(a, v) * (1 - c);
}






inline bool snapAndExpandNormals(EdgeDualNormalGeometry& geometry,
                                 float angleTolerance = 25.0f,
                                 float expandAngle = 5.0f) {
    auto normals1 = geometry.getNormals1();
    auto normals2 = geometry.getNormals2();
    if (normals1.size() != normals2.size()) return false;

    float cosTol = cos(angleTolerance * M_PI / 180.0f);
    float expandRad = expandAngle * M_PI / 180.0f;

    for (size_t i = 0; i < normals1.size(); ++i) {
        Vector3 n1 = normalize(normals1[i]);
        Vector3 n2 = normalize(normals2[i]);

        // Calculate dot product (cosine of angle between normals)
        float dotp = dot(n1, n2);
        
        // Clamp dot product to [-1, 1] to handle numerical errors
        dotp = std::max(-1.0f, std::min(1.0f, dotp));
        

        if (dotp >= cosTol) { // 角度小于阈值
            Vector3 axis = cross(n1, n2);
            if (norm(axis) < 1e-6) continue; // 平行或反平行，跳过

            axis = normalize(axis);

            // 一个转正，一个转负
            normals1[i] = normalize(rotateAroundAxis(n1, axis, -expandRad));
            normals2[i] = normalize(rotateAroundAxis(n2, axis,  expandRad));
        }
    }

    geometry.setNormals1(normals1);
    geometry.setNormals2(normals2);
    return true;
}


} // namespace geometrycentral
