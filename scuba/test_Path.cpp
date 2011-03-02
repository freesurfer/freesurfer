/**
 * @file  test_Path.cpp
 * @brief test Path class
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.10 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <fstream>
#include <sstream>
#include "Path.h"
#include "Scuba-impl.h"
#include "Listener.h"

const char* Progname = "test_Path";

using namespace std;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  cerr << ss.str().c_str() << endl; \
  throw runtime_error( ss.str() ); \
  }


class PathTester {
public:
  void Test();
};

class PathListener : public Listener {
public:
  PathListener() : Listener("PathListener") {
    mPathID = -1;
    mbPathChanged = false;
    mbPathVertexAdded = false;
  }

  virtual void DoListenToMessage ( string isMessage, void* iData ) {
    if ( isMessage == "pathChanged" ) {
      mbPathChanged = true;
      mPathID = *(int*)iData;
    }
    if ( isMessage == "pathVertexAdded" ) {
      mbPathVertexAdded = true;
      mPathID = *(int*)iData;
    }
  }

  void Reset () {
    mPathID = -1;
    mbPathChanged = false;
    mbPathVertexAdded = false;
  }
  bool IsPathChanged ()     {
    return mbPathChanged;
  }
  bool IsPathVertexAdded () {
    return mbPathVertexAdded;
  }
  int  GetPathID ()         {
    return mPathID;
  }
  bool mbPathChanged;
  bool mbPathVertexAdded;
  int  mPathID;
};

void
PathTester::Test () {

  stringstream ssError;

  try {

    Path<float> l;
    PathListener listener;

    Assert( (l.GetID() == 0), "GetID() didn't work." );
    Path<float>& a = Path<float>::FindByID( l.GetID() );
    Assert( (l.GetID() == a.GetID()), "FindByID() didn't work." );

    l.AddListener( listener );

    Point3<float> f;
    Point3<float>& t = f;

    for ( int c = 0; c < 100; c++ ) {
      f.Set( c, c, c );
      l.AddVertex( f );
      t = l.mVertices[c];
      Assert( (t[0] == f[0] && t[1] == f[1] && t[2] == f[2]),
              "Got a point back but it wasn't equal to the one we put in" );
      Assert( (listener.IsPathVertexAdded() &&
               l.GetID() == listener.GetPathID()),
              "Path did not broadcast proper pathVertexAdded msg." );
      listener.Reset();
    }

    Assert( (l.mVertices.size() == 100),
            "Incorrect size" );

    l.MarkEndOfSegment();

    t = l.GetPointAtEndOfLastSegment();
    Assert( (t[0] == 99 && t[1] == 99 && t[2] == 99),
            "Point at end of last segment not correct" );

    for ( int c = 100; c < 200; c++ ) {
      f.Set( c, c, c );
      l.AddVertex( f );
      t = l.mVertices[c];
      Assert( (t[0] == f[0] && t[1] == f[1] && t[2] == f[2]),
              "Got a point back but it wasn't equal to the one we put in" );
    }

    t = l.GetPointAtEndOfLastSegment();
    Assert( (t[0] == 99 && t[1] == 99 && t[2] == 99),
            "Point at end of last segment not correct" );

    Assert( (l.mVertices.size() == 200),
            "Incorrect size" );

    l.ClearLastSegment();

    {
      stringstream ssError;
      ssError << "Incorrect size after ClearLastSegment, was " 
	      << l.mVertices.size();
      Assert( (l.mVertices.size() == 100), ssError.str() );
    }
            

    for ( int c = 0; c < 100; c++ ) {
      t = l.mVertices[c];
      Assert( (t[0] == c && t[1] == c && t[2] == c),
              "Got a point back but it wasn't equal to the one we put in" );
    }

    for ( int c = 100; c < 200; c++ ) {
      f.Set( c, c, c );
      l.AddVertex( f );
      t = l.mVertices[c];
      Assert( (t[0] == f[0] && t[1] == f[1] && t[2] == f[2]),
              "Got a point back but it wasn't equal to the one we put in" );
    }

    Assert( (l.mVertices.size() == 200),
            "Incorrect size" );


    // Test clear.
    l.Clear();
    Assert( (l.mVertices.size() == 0), "Incorrect size after clear" );
    Assert( (listener.IsPathChanged() && l.GetID() == listener.GetPathID()),
            "Path did not broadcast proper pathChanged msg." );
    listener.Reset();

    // Make some verts and move them.
    for ( int c = 0; c < 10; c++ ) {
      f.Set( c, c, c );
      l.AddVertex( f );
    }
    Point3<float> m( 1, -1, 0 );
    l.Move( m );
    Assert( (listener.IsPathChanged() && l.GetID() == listener.GetPathID()),
            "Path did not broadcast proper pathChanged msg." );
    listener.Reset();
    for ( int c = 0; c < 10; c++ ) {
      t = l.mVertices[c];
      Assert( (t[0] == c+1 && t[1] == c-1 && t[2] == c),
              "Vertices not moved correctly." );
    }



    // Test stream r/writing.
    ofstream sw( "/tmp/path.test", ios::out );
    int const cPaths = 10;
    int const cVerts = 100;
    Path<float> pw[cPaths];
    Path<float> pr[cPaths];
    for ( int nPath = 0; nPath < cPaths; nPath++ ) {
      for ( int nVertex = 0; nVertex < cVerts; nVertex++ ) {
        Point3<float> v( nPath, nVertex, nPath*nVertex );
        pw[nPath].AddVertex( v );
      }

      pw[nPath].WriteToStream( sw );
    }
    sw.close();

    ifstream sr( "/tmp/path.test", ios::in );
    for ( int nPath = 0; nPath < cPaths; nPath++ ) {
      pr[nPath].ReadFromStream( sr );
      if ( (pr[nPath].GetNumVertices() != pw[nPath].GetNumVertices()) ) {
        stringstream ssErr;
        ssErr << "Path didn't read correct num verts from stream: "
        << "wrote " << pw[nPath].GetNumVertices()
        << " read " << pr[nPath].GetNumVertices();
        Assert( 0, ssErr.str() );
      }


      for ( int nVertex = 0; nVertex < cVerts; nVertex++ ) {
        Point3<float> vr = pr[nPath].GetVertexAtIndex( nVertex );
        Point3<float> vw = pw[nPath].GetVertexAtIndex( nVertex );
        Assert( (vr[0] == vw[0] && vr[1] == vw[1] && vr[2] == vw[2]),
                "Vertices read from stream didn't match vertices written." );
      }
    }
    sr.close();

  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }
}


int main ( int argc, char** argv ) {

  cerr << "Beginning test" << endl;

  try {

    PathTester tester0;
    for ( int i = 0; i < 1; i++ ) {
      tester0.Test();
    }

  } catch ( runtime_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  } catch (...) {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

