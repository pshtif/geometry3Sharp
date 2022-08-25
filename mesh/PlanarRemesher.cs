using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace g3.mesh
{
    public class PlanarRemesher : MeshRefinerBase
    {
        public PlanarRemesher(DMesh3 m) : base(m)
        {

        }

        /// <summary>
        /// Number of edges that were modified in previous Remesh pass.
        /// If this number gets small relative to edge count, you have probably converged (ish)
        /// </summary>
        public int ModifiedEdgesLastPass = 0;


        protected enum ProcessResult
        {
            Ok_Flipped,
            Ok_Removed,
            Ignored_EdgeIsFine,
            Ignored_EdgeIsFullyConstrained,
            Failed_OpNotSuccessful,
            Failed_NotAnEdge
        };


        public void Remesh()
        {
            // we need to process removal once, then iterate flip and removal until useful
            var removed = MergeTrianglesMultiPass();
            while (true)
            {
                EdgeFlipPass();
                if (ModifiedEdgesLastPass == 0)
                    break;
                removed = MergeTrianglesMultiPass();
                if (removed == 0)
                    break;
            }
            mesh.CleanupUnusedVertices();
        }

        public int MergeTrianglesMultiPass()
        {
            int total = 0;
            do
            {
                MergeTrianglesPass();
                total += ModifiedEdgesLastPass;
            } while (ModifiedEdgesLastPass > 0);
            return total;
        }

        public void MergeTrianglesPass()
        {
            bool done;
            ModifiedEdgesLastPass = 0;
            int cur_eid = start_edges();
            do
            {
                if (mesh.IsEdge(cur_eid))
                {
                    ProcessResult result = ProcessEdgeRemoval(cur_eid);
                    if (result == ProcessResult.Ok_Removed)
                        ModifiedEdgesLastPass++;
                }
                if (Cancelled())        
                    return;
                cur_eid = next_edge(cur_eid, out done);
            } while (done == false);
        }

        double precision { get; set; } = 0.0001;

        private ProcessResult ProcessEdgeRemoval(int edgeID)
        {
            EdgeConstraint constraint =
               (constraints == null) ? EdgeConstraint.Unconstrained : constraints.GetEdgeConstraint(edgeID);
            if (constraint.NoModifications)
                return ProcessResult.Ignored_EdgeIsFullyConstrained;

            // look up verts and tris for this edge
            int a = 0, b = 0, t0 = 0, t1 = 0;
            if (mesh.GetEdge(edgeID, ref a, ref b, ref t0, ref t1) == false)
                return ProcessResult.Failed_NotAnEdge;

            bool bIsBoundaryEdge = (t1 == DMesh3.InvalidID);
            if (bIsBoundaryEdge)
                return ProcessResult.Ignored_EdgeIsFine;

            

            Index2i ov = mesh.GetEdgeOpposingV(edgeID);
            int c = ov.a, d = ov.b;


            var removalCandidates = new List<int>() { a, b }.ToArray();
            var oppositesVector = mesh.GetVertex(d) - mesh.GetVertex(c);
            oppositesVector.Normalize();

            var keep = new List<int>(2);
            foreach (var removalCandidate in removalCandidates)
            {
                var candidateVector = mesh.GetVertex(removalCandidate) - mesh.GetVertex(c);
                var proj = oppositesVector * oppositesVector.Dot(candidateVector);
                if (!proj.EpsilonEqual(candidateVector, precision))
                    keep.Add(removalCandidate);
            }

            if (keep.Count == 1) // only 1 vertex to keep, we can collapse the two triangles to one
            {
                // we can remove the two involved triangles and replace for a single one
                var rem1 = mesh.RemoveTriangle(t0, false, false);
                if (rem1 != MeshResult.Ok)
                    return ProcessResult.Failed_OpNotSuccessful;
                var rem2 = mesh.RemoveTriangle(t1, false, false);
                if (rem2 != MeshResult.Ok)
                    return ProcessResult.Failed_OpNotSuccessful;
                var add = mesh.AppendTriangle(c, d, keep[0]);
                if (add == DMesh3.InvalidID)
                    return ProcessResult.Failed_OpNotSuccessful;
                return ProcessResult.Ok_Removed;
            }
            return ProcessResult.Ignored_EdgeIsFine;
        }

        public void EdgeFlipPass()
        {
            PreventFlip = new List<int>();
            bool done;
            ModifiedEdgesLastPass = 0;
            int cur_eid = start_edges();
            do
            {
                if (mesh.IsEdge(cur_eid))
                {
                    ProcessResult result = ProcessEdgeFlip(cur_eid);
                    if (result == ProcessResult.Ok_Flipped)
                        ModifiedEdgesLastPass++;
                }
                if (Cancelled())        // expensive to check every iter?
                    return;
                cur_eid = next_edge(cur_eid, out done);
            } while (done == false);
        }

        List<int> PreventFlip = new List<int>();

        // most of this boilerplate comes from the Remesher class


        private ProcessResult ProcessEdgeFlip(int edgeID)
        {
            EdgeConstraint constraint =
                (constraints == null) ? EdgeConstraint.Unconstrained : constraints.GetEdgeConstraint(edgeID);
            if (constraint.NoModifications)
                return ProcessResult.Ignored_EdgeIsFullyConstrained;
            if (PreventFlip.Contains(edgeID))
                return ProcessResult.Ignored_EdgeIsFine;
            if (!constraint.CanFlip)
                return ProcessResult.Ignored_EdgeIsFine;


            // look up verts and tris for this edge
            int a = 0, b = 0, t0 = 0, t1 = 0;
            if (mesh.GetEdge(edgeID, ref a, ref b, ref t0, ref t1) == false)
                return ProcessResult.Failed_NotAnEdge;
            bool bIsBoundaryEdge = (t1 == DMesh3.InvalidID);
            if (bIsBoundaryEdge)
                return ProcessResult.Ignored_EdgeIsFine;

            // avoid to flip edge if the shape of t1 and t2 is concave so that 
            // external the profile of the flip stays the same
            //
            Index2i ov = mesh.GetEdgeOpposingV(edgeID);
            int c = ov.a, d = ov.b;
            if (flip_inverts_normals(a, b, c, d, t0))
                return ProcessResult.Ignored_EdgeIsFine;


            // only flip if two trinagles share normal
            // mesh.GetTriangle(t0).
            var tri0 = mesh.GetTriangle(t0);
            Vector3d n0 = MathUtil.Normal(
                mesh.GetVertex(tri0.a),
                mesh.GetVertex(tri0.b),
                mesh.GetVertex(tri0.c)
                );
            var tri1 = mesh.GetTriangle(t1);
            Vector3d n1 = MathUtil.Normal(
                mesh.GetVertex(tri1.a),
                mesh.GetVertex(tri1.b),
                mesh.GetVertex(tri1.c)
                );

            if (!n1.EpsilonEqual(n0, precision))
                return ProcessResult.Ignored_EdgeIsFine;

            var thisFlips = new List<int>();
            thisFlips.AddRange(mesh.GetTriEdges(t0).array);
            thisFlips.AddRange(mesh.GetTriEdges(t1).array);


            DMesh3.EdgeFlipInfo flipInfo;
            MeshResult result = mesh.FlipEdge(edgeID, out flipInfo);
            if (result == MeshResult.Ok)
            {
                PreventFlip.AddRange(thisFlips);
                // DoDebugChecks();
                return ProcessResult.Ok_Flipped;
            }
            return ProcessResult.Failed_OpNotSuccessful;
        }



        // We are using a modulo-index loop to break symmetry/pathological conditions. 
        // For example in a highly tessellated minimal cylinder, if the top/bottom loops have
        // sequential edge IDs, and all edges are < min edge length, then we can easily end
        // up successively collapsing each tiny edge, and eroding away the entire mesh!
        // By using modulo-index loop we jump around and hence this is unlikely to happen.
        const int nPrime = 31337;     // any prime will do...
        int nMaxEdgeID;
        protected virtual int start_edges()
        {
            nMaxEdgeID = mesh.MaxEdgeID;
            return 0;
        }

        protected virtual int next_edge(int cur_eid, out bool bDone)
        {
            int new_eid = (cur_eid + nPrime) % nMaxEdgeID;
            bDone = (new_eid == 0);
            return new_eid;
        }

    }
}
