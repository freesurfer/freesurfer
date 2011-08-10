#ifndef _CDT_H_
#define _CDT_H_
#include <stdio.h>


// aside from small modifications, this code is by Dani Lischinski (http://www.cs.huji.ac.il/~danix/)


/*****************************************************************************
**
** Copyright (C) 1995 by Dani Lischinski 
**
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
**
******************************************************************************/


#define VERBOSE 0
#define NIL(x) NULL
#include <assert.h>
#ifndef MAX
	#define MAX( a, b ) (a) > (b) ? (a) : (b)
#endif
#define MARKED 1


//-------------------------------------------
// FROM "common.h>
//-------------------------------------------


inline int imin(double a, double b, double c)
{
	return (a < b) ? (a < c ? 0 : 2)
				   : (b < c ? 1 : 2);
}


//-------------------------------------------
// FROM "llist.h>
//-------------------------------------------


typedef void* Ldata;

class Llist;

class LlistNode {
	friend class Llist;
private:
	Ldata data;
	LlistNode* next;
	LlistNode() 						{ data = 0; next = this; }
	LlistNode(Ldata d, LlistNode* n)	{ data = d; next = n; }
};

typedef LlistNode *LlistPos;

class Llist {
private:
	int len;
	LlistNode *node;
public:
	Llist()						{ node = new LlistNode; len = 0; }
	~Llist();
	bool empty() const		{ return (node == node->next); }
	void insert(LlistPos, Ldata);
	void remove(LlistPos);
	LlistPos first() const			{ return node; }
	bool isEnd(LlistPos p) const	{ return (p->next == node); }
	LlistPos next(LlistPos p) const	{ return p->next; }
	Ldata retrieve(LlistPos p) const	{ return p->next->data; }
	void store(LlistPos p, Ldata d)	{ p->next->data = d; }
	int length() const			{ return len; }
};

inline void Llist::insert(LlistPos p, Ldata d)
{
	p->next = new LlistNode(d, p->next);
	len++;
}

inline void Llist::remove(LlistPos p)
{
	LlistNode *t = p->next;
	p->next = p->next->next;
	delete t;
	len--;
}

inline Llist::~Llist()
{
	while (!empty()) remove(first());
	delete node;
}


//-------------------------------------------
// FROM "geom2d.h>
//-------------------------------------------


#include <math.h>


#define EPS 1e-6
#define X 0
#define Y 1
#define Z 2

typedef double Real;
typedef bool Boolean;

class Vector2d {
public:
	Real x, y;
	void *tag;
	Vector2d()				{ x = 0; y = 0; tag = 0; }
	Vector2d(Real a, Real b)		{ x = a; y = b; tag = 0; }
	Vector2d(Real a, Real b, void *t)	{ x = a; y = b; tag = t; }
	Vector2d(const Vector2d &v)		{ x = v.x; y = v.y; tag = v.tag; }
	Real& operator[](int i)			{ return ((Real *)this)[i]; }
	const Real& operator[](int i) const	{ return ((Real *)this)[i]; }
	Real norm() const;
	void normalize();
	bool operator==(const Vector2d&) const;
	Vector2d operator+(const Vector2d&) const;
	Vector2d operator-(const Vector2d&) const;
	Real     operator|(const Vector2d&) const;
	friend Vector2d operator*(Real, const Vector2d&);
	friend Vector2d operator/(const Vector2d&, Real);
//	friend istream& operator>>(istream&, Vector2d&);
//	friend ostream& operator<<(ostream&, const Vector2d&);
};

typedef Vector2d Point2d;

class Line {
public:
	Line()	{}
	Line(const Point2d&, const Point2d&);
	void set(const Point2d&, const Point2d&);
	Real eval(const Point2d&) const;
	int classify(const Point2d&) const;
	Point2d intersect(const Point2d&, const Point2d&) const;
private:
	Real a, b, c;
};

// Vector2d:

inline Real Vector2d::norm() const
{
	return sqrt(x * x + y * y);
}

inline void Vector2d::normalize()
{
	Real len;

	if ((len = sqrt(x * x + y * y)) == 0.0)
		printf( "Vector2d::normalize: Division by 0\n" );
	else {
		x /= len;
		y /= len;
	}
}

inline Vector2d Vector2d::operator+(const Vector2d& v) const
{
	return Vector2d(x + v.x, y + v.y);
}

inline Vector2d Vector2d::operator-(const Vector2d& v) const
{
	return Vector2d(x - v.x, y - v.y);
}

inline bool Vector2d::operator==(const Vector2d& v) const
{
	return (*this - v).norm() <= EPS;
}

inline Real Vector2d::operator|(const Vector2d& v) const
{
	return x * v.x + y * v.y;
}

inline Vector2d operator*(Real c, const Vector2d& v)
{
	return Vector2d(c * v.x, c * v.y);
}

inline Vector2d operator/(const Vector2d& v, Real c)
{
	return Vector2d(v.x / c, v.y / c);
}

// Line:

inline Line::Line(const Point2d& p, const Point2d& q)
// Computes the normalized line equation through the
// points p and q.
{
	Vector2d t = q - p;
	Real len = t.norm();

	a =   t.y / len;
	b = - t.x / len;
	// c = -(a*p.x + b*p.y);

	// less efficient, but more robust -- seth.
	c = -0.5 * ((a*p.x + b*p.y) + (a*q.x + b*q.y));
}

inline void Line::set(const Point2d& p, const Point2d& q)
{
	*this = Line(p, q);
}

inline Real Line::eval(const Point2d& p) const
// Plugs point p into the line equation.
{
	return (a * p.x + b* p.y + c);
}

inline Point2d Line::intersect(const Point2d& p1, const Point2d& p2) const
// Returns the intersection of the line with the segment (p1,p2)
{
        // assumes that segment (p1,p2) crosses the line
        Vector2d d = p2 - p1;
        Real t = - eval(p1) / (a*d[X] + b*d[Y]);
        return (p1 + t*d);
}

inline int Line::classify(const Point2d& p) const
// Returns -1, 0, or 1, if p is to the left of, on,
// or right of the line, respectively.
{
	Real d = eval(p);
	return (d < -EPS) ? -1 : (d > EPS ? 1 : 0);
}

inline Boolean operator==(const Point2d& point, const Line& line)
// Returns TRUE if point is on the line (actually, on the EPS-slab
// around the line).
{
	Real tmp = line.eval(point);
	return(fabs(tmp) <= EPS);
}

inline Boolean operator<(const Point2d& point, const Line& line)
// Returns TRUE if point is to the left of the line (left to
// the EPS-slab around the line).
{
	return (line.eval(point) < -EPS);
}


//-------------------------------------------
// FROM "quadedge.h>
//-------------------------------------------


#define TRUE true
#define FALSE false

class QuadEdge;

class Edge {
	friend class QuadEdge;
	friend void Splice(Edge*, Edge*);
  private:
	int num;
	Edge *next;
	Point2d *data;
#if MARKED
	short _lmark;
	short _rmark;
#endif
  public:
	Edge()	{ data = 0; }
	Edge* Rot();
	Edge* invRot();
	Edge* Sym();
	Edge* Onext();
	Edge* Oprev();
	Edge* Dnext();
	Edge* Dprev();
	Edge* Lnext();	
	Edge* Lprev();
	Edge* Rnext();	
	Edge* Rprev();
	Point2d* Org();
	Point2d* Dest();
	const Point2d& Org2d() const;
	const Point2d& Dest2d() const;
	void  EndPoints(Point2d*, Point2d*);
	Boolean isConstrained();
	void Constrain();
	QuadEdge* QEdge() const				{ return (QuadEdge *)(this - num); }
#if MARKED
	short& lmark( void ) { return _lmark; }
	short& rmark( void ) { return _rmark; }
#endif
};

class QuadEdge {
  private:
	Edge e[4];
	Boolean c;
#if MARKED
	short _lmark;
	short _rmark;
#endif

  public:
	QuadEdge(Boolean);
	Edge *edges()			{ return e; }
	Boolean isConstrained()	{ return c; }
	void Constrain() 		{ c = TRUE; }
	void Draw();
	~QuadEdge();
#if MARKED
	short& lmark( void ) { return _lmark; }
	short& rmark( void ) { return _rmark; }
#endif
};

class Mesh {
  private:
	Edge *startingEdge;
//	Llist edges;
	void DeleteEdge(Edge*);
	Edge *Connect(Edge*, Edge*);
	Edge *Locate(const Point2d&);
	Edge *BrutalLocate(const Point2d&);
	void SplitEdge(Edge*, const Point2d&);
	void Triangulate(Edge*);
  public:
	Llist edges;
	Mesh(const Point2d&, const Point2d&, const Point2d&);
	Mesh(const Point2d&, const Point2d&, const Point2d&, const Point2d&);
	Edge *MakeEdge(Boolean);
	Edge *MakeEdge(Point2d*, Point2d*, Boolean);
	Edge *InsertSite(const Point2d&, Real dist = EPS);
	void InsertEdge(const Point2d&, const Point2d&);
	int  numEdges() const	{ return edges.length(); }
	void Draw() const;
	void Apply(void (*f)(int isconstrained, double x1, double y1, double x2, double y2)) const;
#if MARKED
	void ApplyTriangles(void (*f)( Point2d p1, Point2d p2, Point2d p3 )) const;
#endif
	void ApplyTaggedMesh(void (*f)(int isconstrained, void *tag1, void *tag2)) const;
	void ApplyTaggedKey(void *key1, void *key2,
			    int (*match)(void *key, void *tag),
			    void (*f)(int isconstrained, void *key1, void *key2, void *tag1, void *tag2)) const;
	// true if no edge in mesh matches given tags
	int TaggedEdgeUnmatched(void *key1, void *key2, int (*match)(void *key, void *tag)) const;

	// destructor
	~Mesh();
};

inline QuadEdge::QuadEdge(Boolean constrained = FALSE)
{
	e[0].num = 0, e[1].num = 1, e[2].num = 2, e[3].num = 3;
	e[0].next = &(e[0]); e[1].next = &(e[3]);
	e[2].next = &(e[2]); e[3].next = &(e[1]);
	c = constrained;
}

/************************* Edge Algebra *************************************/

inline Edge* Edge::Rot()
// Return the dual of the current edge, directed from its right to its left. 
{
	return (num < 3) ? this + 1 : this - 3;
}

inline Edge* Edge::invRot()
// Return the dual of the current edge, directed from its left to its right. 
{
	return (num > 0) ? this - 1 : this + 3;
}

inline Edge* Edge::Sym()
// Return the edge from the destination to the origin of the current edge.
{
	return (num < 2) ? this + 2 : this - 2;
}
	
inline Edge* Edge::Onext()
// Return the next ccw edge around (from) the origin of the current edge.
{
	return next;
}

inline Edge* Edge::Oprev()
// Return the next cw edge around (from) the origin of the current edge.
{
	return Rot()->Onext()->Rot();
}

inline Edge* Edge::Dnext()
// Return the next ccw edge around (into) the destination of the current edge.
{
	return Sym()->Onext()->Sym();
}

inline Edge* Edge::Dprev()
// Return the next cw edge around (into) the destination of the current edge.
{
	return invRot()->Onext()->invRot();
}

inline Edge* Edge::Lnext()
// Return the ccw edge around the left face following the current edge.
{
	return invRot()->Onext()->Rot();
}

inline Edge* Edge::Lprev()
// Return the ccw edge around the left face before the current edge.
{
	return Onext()->Sym();
}

inline Edge* Edge::Rnext()
// Return the edge around the right face ccw following the current edge.
{
	return Rot()->Onext()->invRot();
}

inline Edge* Edge::Rprev()
// Return the edge around the right face ccw before the current edge.
{
	return Sym()->Onext();
}

/************** Access to non-topological info ******************************/

inline Point2d* Edge::Org()
{
	return data;
}

inline Point2d* Edge::Dest()
{
	return Sym()->data;
}

inline const Point2d& Edge::Org2d() const
{
	return *data;
}

inline const Point2d& Edge::Dest2d() const
{
	return (num < 2) ? *((this + 2)->data) : *((this - 2)->data);
}

inline void Edge::EndPoints(Point2d* orig, Point2d* de)
{
	data = orig;
	Sym()->data = de;
}

inline Boolean Edge::isConstrained()
{
	return QEdge()->isConstrained();
}

inline void Edge::Constrain()
{
	QEdge()->Constrain();
}


#endif // _CDT_H_
