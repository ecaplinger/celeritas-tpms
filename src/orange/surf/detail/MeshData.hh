#pragma once

#include "corecel/cont/Array.hh"
#include "corecel/data/Collection.hh"

namespace celeritas
{

struct Triangle
{
    Array<int, 3> verts;

    Triangle(Array<int, 3> verts)
        : verts(verts)
    {
    }

    int operator[](size_type idx) const { return verts[idx]; }
};

template<Ownership W, MemSpace M>
struct MeshData
{
    Collection<Real3, W, M> global_verts;
    Collection<Triangle, W, M> triangles;

    Array<Real3, 3> triangle_at(int idx)
    {
        const Triangle& tri = triangles[ItemId<Triangle>(idx)];
        return {
            global_verts[ItemId<Real3>(tri[0])],
            global_verts[ItemId<Real3>(tri[1])],
            global_verts[ItemId<Real3>(tri[2])],
        };
    }
};

template<MemSpace M>
struct DirectionalMeshData
{
    Collection<Real3, Ownership::value, M> global_verts;
    Collection<Real3, Ownership::value, M> global_normals;      // TODO: evaluate whether this is necessary
    Collection<Triangle, Ownership::value, M> triangles;
};

}
