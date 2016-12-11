﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace g3
{
    public class OBJWriter : IMeshWriter
    {
        public IOWriteResult Write(BinaryWriter writer, List<IMesh> vMeshes, WriteOptions options)
        {
            // [RMS] not supported
            throw new NotImplementedException();
        }



        public IOWriteResult Write(TextWriter writer, List<IMesh> vMeshes, WriteOptions options)
        {
            int nAccumCountV = 1;       // OBJ indices always start at 1

            for (int mi = 0; mi < vMeshes.Count; ++mi) {

                IMesh mesh = vMeshes[mi];
                bool bVtxColors = options.bPerVertexColors && mesh.HasVertexColors;
                bool bNormals = options.bPerVertexNormals && mesh.HasVertexNormals;

				int[] mapV = new int[mesh.MaxVertexID];

                foreach ( int vi in mesh.VertexIndices() ) { 
					mapV[vi] = nAccumCountV++;
                    Vector3d v = mesh.GetVertex(vi);
                    if ( bVtxColors ) {
                        Vector3d c = mesh.GetVertexColor(vi);
                        writer.WriteLine("v {0} {1} {2} {3:F8} {4:F8} {5:F8}", v[0], v[1], v[2], c[0],c[1],c[2]);
                    } else {
                        writer.WriteLine("v {0} {1} {2}", v[0], v[1], v[2]);
                    }

                    if ( options.bPerVertexNormals && mesh.HasVertexNormals ) {
                        Vector3d n = mesh.GetVertexNormal(vi);
                        writer.WriteLine("vn {0:F10} {1:F10} {2:F10}", n[0], n[1], n[2]);
                    }
                }

                foreach (int ti in mesh.TriangleIndices() ) { 
                    Vector3i t = mesh.GetTriangle(ti);
					t[0] = mapV[t[0]];
					t[1] = mapV[t[1]];
					t[2] = mapV[t[2]];
                    
                    if ( bNormals ) {
                        writer.WriteLine("f {0}//{0} {1}//{1} {2}//{2}", t[0], t[1], t[2]);
                    } else {
                        writer.WriteLine("f {0} {1} {2}", t[0], t[1], t[2]);
                    }

                }

            }


            return new IOWriteResult(WriteResult.Ok, "");
        }


    }
}
