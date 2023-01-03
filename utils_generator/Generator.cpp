// Compile with     g++ -std=c++11 Generator.cpp
// Run with         rm -rf tmp ; mkdir tmp ; ./a.out ./tmp/
// Merge with       bcompare ./tmp ../include &
//                  cp tmp/* ../include
//
// Abstractly, we have
//			A Surface,			with various properties
//			A set of vertices,	with various properties
//			A set of faces,		with various properties
// We have several different representations of this
//			The traditional MRIS struct, containing pointer to vectors of VERTEX_TOPOLOGY, VERTEX and FACE
//				This is fully implemented in the Freesurfer code
//			The MRIS_MP struct, which contains pointers to vectors of some of the vertices and faces properties, but not all
//			and which can get some, but not all, of the other properties from an underlying MRIS struct
//				This is fully implemented in the Freesurfer code, being initialized from an MRIS struct
//			The MRISPV struct, which which contains pointers to vectors of all of the vertices and faces properties
//				This is not yet implemented in the code
//
// For each of these representations, we want to generate a cascade of Surface,Face,Vertex accessing classes that provide
// only limited access, so we can try to restrict phases of the Surface construction and manipulation to only parts of the 
// whole data structure to force the data structures to be internally consistent by restricting who can modify it, and to 
// try to establish a dependency order onto the properties.
//
// To achieve this, the following code is composed of several sections
//
//		Section 1: Utilities					Useful classes and functions that are mostly independent of the goal
//		Section 2: Abstract	representation		Classes etc. that describe the abstraction representation
//												These are used to create a Representations - the complete list of properties etc.
//		Section 3: Generator utilities			Useful classes etc. that all the source generators share
//		Section 4: Representation Generators	One for each of the representations described above
//		Section 5: Accessor Generators			One for each of the representations described above
//		Section 6: main							Invoke some of the generators
//		Section 7: Representations						Define the properties, what each representation supports, what each phase can access


// Section 1: Utilities		Useful classes and functions that are mostly independent of the goal
//
#include <iostream>
#include <iomanip>
#include <fstream>

#include <algorithm>
#include <list>
#include <locale>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <assert.h>

static std::string uppercase(std::string const & s_init) {
	using namespace std;
	string s = s_init;
	transform(s.begin(), s.end(), s.begin(), ::toupper);
	return s;
}

static std::ostream& indent(std::ostream & os, size_t depth) {
	for (size_t i = 0; i < depth; i++) os << "    ";
	return os;
}

void bpt() {}

namespace ColumnedOutput {

	using namespace std;

	enum Justify {Just_left, Just_right, Just_none};

	struct T {
		
		struct Cell {
			string s;
			bool ignoreWidth;
			Justify justify;
			Cell() : ignoreWidth(false), justify(Just_none) {}
		};

		T() : col(0), newRowPending(true) {}

		Cell & cell() {
			if (newRowPending) newRow();
			auto & r = grid.back();
			while (col >= r.size()) r.push_back(Cell());
			return r[col];
		}

		T& operator <<(std::string const & s) {
			auto & c = cell();
			c.s += s;
			if (!c.ignoreWidth) {
				while (colSizes.size() <= col) colSizes.push_back(0);
				colSizes[col] = max(colSizes[col], c.s.size());
			}
			return *this;
		}

		struct IgnoreWidth{};
		T& operator <<(IgnoreWidth const & s) {
			cell().ignoreWidth = true;
			return *this;
		}

		T& operator <<(Justify to) {
			cell().justify = to;
			return *this;
		}

		struct EndC {};
		T& operator <<(EndC const &) {
			col++;
			return *this;
		}

		struct EndR {};
		T& operator <<(EndR const &) {
			newRowPending = true;
			col = 0;
			return *this;
		}

		void append(ostream & os, size_t depth) {
			for (auto const & r : grid) {
				indent(os,depth);
				auto sep = "";
                
                auto r_size = r.size();
                while (r_size && (r[r_size-1].s.find_first_not_of(" \t") == string::npos)) r_size--;   // avoid trailing blanks
                                
				for (size_t ci = 0; ci < r_size; ci++) {
					os << sep;
                    if (ci+1 < r_size || r[ci].justify == Just_right) {
					    if (ci < colSizes.size() && colSizes[ci] > 0)
						    os << setw(colSizes[ci]);
                    }
					switch (r[ci].justify) {
					case Just_left:  os << left;  break;
					case Just_right: os << right; break;
					}
					os << r[ci].s; sep = " ";
				}
				os << endl;
				os << left;
			}
		}

	private:
		bool newRowPending;
		size_t col;
		vector< vector<Cell> > grid;
		vector<size_t> colSizes;
		void newRow() {
			grid.push_back(vector<Cell>());
			newRowPending = false;
		}
	};

	T::IgnoreWidth ignoreWidth;
	T::EndC endC;
	T::EndR endR;
};

static std::string createFilesWithin;


// Section 2: Abstract	representation		Classes etc. that describe the abstraction representation
//
namespace Abstract_Representation {
	using namespace std;
	using namespace ColumnedOutput;

#define PHASES \
	ELT(ExistenceM 				, Existence					) SEP \
	ELT(Existence				, end						) SEP \
	ELT(TopologyM				, Topology					) SEP \
	ELT(Topology				, Existence					) SEP \
	ELT(XYZPositionM			, XYZPosition				) SEP \
	ELT(XYZPosition				, Topology					) SEP \
	ELT(XYZPositionConsequencesM, XYZPositionConsequences	) SEP \
	ELT(XYZPositionConsequences	, XYZPosition				) SEP \
	ELT(DistortM				, Distort					) SEP \
	ELT(Distort					, XYZPositionConsequences	) SEP \
	ELT(AnalysisM				, Analysis					) SEP \
	ELT(Analysis				, Distort					) SEP \
	ELT(AllM					, end                       )
	// end of macro

	namespace Phase {
		enum T {
#define SEP 
#define ELT(N,ENHANCES) N,
			PHASES
#undef ELT
#undef SEP
			end
		};
		const char* spelling[] = {
#define SEP ,
#define ELT(N,ENHANCES) #N
			PHASES
#undef ELT
#undef SEP
		};
		string namespaceName(T p) { return spelling[p]; }
		const Phase::T enhances[Phase::end] = {
#define SEP ,
#define ELT(N,ENHANCES) ENHANCES
			PHASES
#undef ELT
#undef SEP
		};
	};

	typedef string Id;

	struct Prop;
	struct Fmbr;
	struct Type;

	enum CommentNature { ComGeneral, ComList, ComListSublist };

	struct PropModifiable {
		CommentNature commentNature;
		bool  nohash;
		string which;
		Prop* repeatedSize;	// type must be a PointerToRepeatedType type, this is the element that gives the number of repeats
		PropModifiable() : commentNature(ComGeneral), nohash(false), repeatedSize(nullptr) {}
	};
	struct Prop : public PropModifiable {
		string      const accessorClassId;
		Type*		const type;
		Id			const id;

        string key() const {
            return accessorClassId+"."+id;
        }

		string      const comment;
		Phase::T    const firstReadingPhase;
		Phase::T	const firstWritingPhase;
		Phase::T	const lastWritingPhase;
		Prop(string const & accessorClassId, Type* t, Id i, Phase::T wb, string const & com = ""                          ) : accessorClassId(accessorClassId), type(t), id(i), firstReadingPhase(wb), firstWritingPhase(wb), lastWritingPhase(wb), comment(com) {}
		Prop(string const & accessorClassId, Type* t, Id i, Phase::T wb, Phase::T we, string const & com = ""             ) : accessorClassId(accessorClassId), type(t), id(i), firstReadingPhase(wb), firstWritingPhase(wb), lastWritingPhase(we), comment(com) {}
		Prop(string const & accessorClassId, Type* t, Id i, Phase::T rb, Phase::T wb, Phase::T we, string const & com = "") : accessorClassId(accessorClassId), type(t), id(i), firstReadingPhase(rb), firstWritingPhase(wb), lastWritingPhase(we), comment(com) {}

		bool isGeneralComment() const { return !type && commentNature == ComGeneral; }

		Prop* setCommentNature(CommentNature to) { commentNature = to;     return this; }
		Prop* setNoHash() { nohash = true;     return this; }
		Prop* setPRSize(Prop* to) { repeatedSize = to; return this; }
		Prop* setWhich(string const & to) { which = to; return this; }
	};

	Phase::T phaseRbegin = Phase::end;
	Phase::T phaseWbegin = Phase::end;
	Phase::T phaseWend   = Phase::end;

	struct AtomicType; struct PointerType; struct ArrayOfPointerType; struct PointerToRepeatedAtomicType;
	struct Type {
		Id			const id;
		Type(Id i) : id(i) {}
		virtual AtomicType*  toAtomicType() { return nullptr; }
		virtual PointerType* toPointerType() { return nullptr; }
		virtual ArrayOfPointerType*  toArrayOfPointerType() { return nullptr; }
		virtual PointerToRepeatedAtomicType* toPointerToRepeatedType() { return nullptr; }
	};

	struct AtomicType : public Type {
		virtual AtomicType*  toAtomicType() { return this; }
		AtomicType(Id i) : Type(i) {}
	};

	struct PointerType : public Type {
		virtual PointerType*  toPointerType() { return this; }
		Type* const target;
		PointerType(Id i, Type* t) : Type(i), target(t) {}
	};

	struct ArrayOfPointerType : public Type {
		virtual ArrayOfPointerType*  toArrayOfPointerType() { return this; }
		ArrayOfPointerType(Id i) : Type(i) {}
	};

	struct PointerToRepeatedAtomicType : public PointerType {
		virtual PointerToRepeatedAtomicType* toPointerToRepeatedType() { return this; }
		PointerToRepeatedAtomicType(Id i, AtomicType* at) : PointerType(i, at) {}
	};

	struct Representation;

	int numericPrecision(string const & t) {
		if (t == "char" || t == "uchar") return 0;
		if (t == "short") return 1;
		if (t == "int") return 2;
		if (t == "float") return 3;
		if (t == "size_t") return 4;
		return -1;
	}
	bool cvtNeededToFrom(string const & to, string const & from) {
		if (from == to) return false;
		auto from_prec = numericPrecision(from);
		auto to_prec = numericPrecision(to);
		if (from_prec < 0 || to_prec < 0) return true;
		return from_prec > to_prec;
	}

	struct How {
		virtual ~How() {}
		virtual bool   isImplDetail ()					const = 0;
		virtual string secondaryName()					const = 0;
		virtual string reprNoArrow	()					const = 0;
		virtual string memberType	(Prop const & prop) const = 0;
		virtual string memberId		(Prop const & prop) const = 0;
		virtual string storeType	(Prop const & prop) const = 0;
		virtual string storeId		(Prop const & prop) const = 0;
		virtual string retType		(Prop const & prop) const = 0;
		virtual string getId		(Prop const & prop) const = 0;
		virtual string getArgs		(Prop const & prop) const = 0;
		virtual string setId		(Prop const & prop) const = 0;
		virtual string indexArg		(Prop const & prop) const = 0;
		virtual string setToArg		(Prop const & prop) const = 0;
		virtual string getSetExpr	(string const & reprNoArrow, Prop const & prop, bool write) const = 0;
	};

	struct HowRedirect : public How {
		How& to;
		HowRedirect(How * to) : to(*to) {}
		HowRedirect(How & to) : to( to) {}
		virtual bool   isImplDetail ()								const { return to.isImplDetail (); }
		virtual string secondaryName()								const { return to.secondaryName(); }
		virtual string reprNoArrow	()								const { return to.reprNoArrow  (); }
		virtual string memberType	(Prop const & prop)				const { return to.memberType(prop); }
		virtual string memberId		(Prop const & prop)				const { return to.memberId  (prop); }
		virtual string storeType	(Prop const & prop)				const { return to.storeType (prop); }
		virtual string storeId		(Prop const & prop)				const { return to.storeId	(prop); }
		virtual string retType		(Prop const & prop)				const { return to.retType	(prop); }
		virtual string getId		(Prop const & prop)				const { return to.getId	    (prop); }
		virtual string getArgs		(Prop const & prop)				const { return to.getArgs	(prop); }
		virtual string setId		(Prop const & prop)				const { return to.setId	    (prop); }
		virtual string indexArg		(Prop const & prop)				const { return to.indexArg  (prop); }
		virtual string setToArg		(Prop const & prop)				const { return to.setToArg  (prop); }
		virtual string getSetExpr	(string const & reprNoArrow, Prop const & prop, bool write) const { return to.getSetExpr(reprNoArrow, prop, write); }
	};

	struct HowDirect : public How {
		virtual bool   isImplDetail()  const { return m_isImplDetail; }
		virtual string secondaryName() const { return m_secondaryName; }
		bool	m_isImplDetail;
		string	m_secondaryName;		// if not "", then one FACE, VERTEX, etc. showing this property has been placed into some composite child vector
										// maybe this should be the prop of the parent's child vector?

		HowDirect() : m_isImplDetail(false) {}
		virtual ~HowDirect() {}

		virtual string reprNoArrow() const {
			return "repr";
		}

		virtual string memberType(Prop const & prop) const {	// the type of the field in the classId
			return prop.type->id;
		}

		virtual string memberId(Prop const & prop) const {		// the id of field in the classId
			return prop.id;										// it might need to be indexed by [idx]
		}

		virtual string storeType(Prop const & prop) const {
			return memberType(prop);
		}
		virtual string storeId(Prop const & prop) const {
			return memberId(prop);
		}
		virtual string retType(Prop const & prop) const {
			return prop.type->id;
		}
		virtual string getId(Prop const & prop) const {
			return prop.id;
		}
		virtual string getArgs(Prop const & prop) const {
			return "";
		}
		virtual string setId(Prop const & prop) const {
			return "set_" + prop.id;
		}
		virtual string indexArg(Prop const & prop) const {
			if (!prop.type->toPointerToRepeatedType()) return "";
			return "size_t i";
		}
		virtual string setToArg(Prop const & prop) const {
			string s = retType(prop) + " to";
			return s;
		}

		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			// generate the string needed to fetch or store the value as a retType
			//	convert between the retType and the storeType 
			//	the input value for a store is in "retType const & to"
			//
			auto st = storeType(prop);
			auto rt = retType(prop);

			string storeName = reprNoArrow + "->" + storeId(prop);
			if (indexArg(prop).size()) {
				storeName += "[i]";
				if (auto pt = prop.type->toPointerToRepeatedType()) {
					st = pt->target->id;
				}
			}
			string fetchName = (write ? "to" : storeName);

			if (rt == "Vertex" || rt == "Face") {
				if (!write) {
					fetchName = rt + "(repr," + storeName + ")";
					st = rt;
				}
			}

			if (cvtNeededToFrom(write?st:rt,write?rt:st)) {
				if   (!write) fetchName = rt + "("+storeName+")";
				else          fetchName = st + "(to)";
			}
			return write ? (storeName + " = " + fetchName) : (fetchName);
		}

	};

	struct Representation {
		struct PropHow {
			Prop*     const prop;
			How*  const how;
			PropHow() : prop(nullptr),   how(nullptr) {}
			PropHow(Prop* prop, How* how) : prop(prop), how(how) {}
		};
		std::vector<PropHow> implements;
		string const reprClassName;
		const char* const rootHeaderFile;
		Representation(const char* rootHeaderFile, string const & surfaceClassName) : rootHeaderFile(rootHeaderFile), reprClassName(surfaceClassName) {}
	};

	template <typename Callable0, typename Callable1, typename Callable2, typename Callable3>
	void walkClasses(
		Representation & representation,
		Callable0 propHowToClassId,
		Callable1 initClass,
		Callable2 nextMember,
		Callable3 finiClass) {

		set<string> classNames;
		for (auto & propHow : representation.implements) {
			if (!propHow.how) continue;
			classNames.insert(propHowToClassId(propHow));
		}

		auto walkClass = [&](bool primaryClass, string classId, string key)
		{
			auto keyIter = classNames.find(key);
			if (keyIter == classNames.end()) return;
			classNames.erase(keyIter);

			if (!classId.size()) classId = key;
			initClass(classId);

			for (int write = 0; write < 2; write++) {
				for (auto & propHow : representation.implements) {
					auto propClassId = propHowToClassId(propHow);
					if (propClassId != key) continue;
					auto & prop = *propHow.prop;
					auto   how  =  propHow.how;
					nextMember(prop, how, write == 1);
				}
			}

			finiClass(classId);
		};

		// the order here matches old code
		
		walkClass(false, "", "face_type_");
		walkClass(false, "", "VERTEX_TOPOLOGY");
		walkClass(false, "", "vertex_type_");

		walkClass(false, "", "Face");
		walkClass(false, "", "Vertex");

		walkClass(true,  representation.reprClassName, "");

		while (!classNames.empty()) {
			auto className = *classNames.begin();
			walkClass(false, className, className);
		}
	}
}


// Section 3: Generator utilities			Useful classes etc. that all the source generators share.
//
namespace Generator_Utilities {
	using namespace std;
	using namespace ColumnedOutput;
	using namespace Abstract_Representation;

	struct Generate_Base {
		Representation & representation;
		
		Generate_Base(
			ostream & init_os,
			Representation & representation) : os(), representation(representation), depth(0) {
			os.push_back(&init_os);
			tos() << endl;
            indent() << "#pragma once"                               << endl;
			indent() << "// GENERATED SOURCE - DO NOT DIRECTLY EDIT" << endl;
			indent() << "// " << endl;
			indent() << "// =======================================" << endl;
            indent() << "#include \"mrisurf_aaa.h\""                 << endl;
		}

		// A stack of output files so we can spread the generated code
		// across several files for readability
		//
		ostream & tos() { return *os.back(); }
		void os_push(ostream * new_tos) { os.push_back(new_tos); }
		void os_pop() { os.pop_back(); }

		// Support indented output to the top output file
		//
		size_t depth;
		ostream& indent() {
			return ::indent(tos(), depth);
		}

	private:
		list<ostream *> os;	// stack of output files
	};
}


// Section 4: Representation Generators			One for each of the representations described above		
//
namespace Representation_Generators {
	using namespace std;
	using namespace ColumnedOutput;
	using namespace Abstract_Representation;
	using namespace Generator_Utilities;

	struct Generate_Representation : public Generate_Base {
		Generate_Representation(
			ostream & os,
			Representation & representation) : Generate_Base(os, representation) {}
		virtual void generateRepresentationClasses() = 0;
	};

	struct Generate_usingHow : public Generate_Representation {

		Generate_usingHow(
			ostream & os,
			Representation & representation) : Generate_Representation(os, representation)
		{
			if (representation.reprClassName == "MRIS") indent() << "#define SEPARATE_VERTEX_TOPOLOGY" << endl;
			generateRepresentationClasses();
			generateMacros();
		}

		void generateRepresentationClasses()
		{
			ColumnedOutput::T* colsPtr;

			walkClasses(representation,
				[&](Representation::PropHow const & propHow) {
					return propHow.how ? propHow.how->secondaryName() : representation.reprClassName;
				},
				[&](string const & classId) {
					indent() << "struct "
						<< classId
						<< " {" << endl;
					depth++;
					colsPtr = new ColumnedOutput::T;
				},
				[&](Prop& prop, How* how, bool write) {
					if (write || !how) return;
					if (!prop.type && prop.commentNature != ComGeneral) return;

					auto & cols = *colsPtr;

					const char * comment = "//";
					if (prop.type) {
						auto memberType = how->memberType(prop);
						cols << memberType << endC << Just_left << how->memberId(prop) << endC << ";";
						comment = "  //";
						if (prop.repeatedSize) {
							cols << ignoreWidth << comment << " size() is " << prop.repeatedSize->id;
							comment = ".  ";
						}
					}
					if (prop.comment.size()) cols << ignoreWidth << comment << "  " << prop.comment;
					cols << endR;
				},
				[&](string const & classId) {
					colsPtr->append(tos(), depth);
					delete colsPtr; colsPtr = nullptr;
					depth--;
					indent() << "};		// " << classId << endl << endl;
				}
			);
		}

		void generateMacros() {
			if (representation.reprClassName != "MRIS") return;

			bool eltEmitted = false;
			bool endlPending = false;
			bool endMPending = false;
			auto emitEndOfMacroIfNeeded = [&] {
				if (!endMPending) return;
				depth--;
				eltEmitted = false;
				if (endlPending) tos() << " \\" << endl; endlPending = false;
				if (endMPending) indent() << "// end of macro" << endl << endl; endMPending = false;
			};

			string secondaryName;

			for (auto & propHow : representation.implements) {
				auto & prop = *propHow.prop;
				auto & how  = *propHow.how;
				if (secondaryName != how.secondaryName()) {
					secondaryName = how.secondaryName();
					emitEndOfMacroIfNeeded();
				}
				if (prop.isGeneralComment()) continue;

				if (eltEmitted)  { tos() << " SEP";		eltEmitted = false; }
				if (endlPending) { tos() << " \\" << endl; endlPending = false; }

				if (!prop.type) {
					switch (prop.commentNature) {
					case ComList:
						emitEndOfMacroIfNeeded();
						indent() << "#define " << prop.comment << " \\" << endl; eltEmitted = false; endlPending = false; endMPending = true;
						depth++;
						break;
					case ComListSublist:
						indent() << prop.comment;  eltEmitted = true; endlPending = true;
						break;
					}
				}
				else if (endMPending) {
					const char* macro = (prop.nohash) ? "ELTX" : "ELTT";
					auto t = prop.type;
					if (auto tp = t->toPointerType()) { macro = "ELTP"; t = tp->target; }
					indent() << macro << "(" << t->id << "," << prop.id << ") "; eltEmitted = true; endlPending = true;
				}
			}

			emitEndOfMacroIfNeeded();
		}

	};
}


// Section 5: Accessor Generators				One for each of the representations described above		
//
namespace Accessor_Generators {
	using namespace std;
	using namespace ColumnedOutput;
	using namespace Abstract_Representation;
	using namespace Generator_Utilities;
	using namespace Representation_Generators;

	struct Generate_abstract_spec_Base;
    struct Generate_abstract_inco_Base;

	struct Generate_Surface_inco;
	struct Generate_Surface_spec;
	struct Generate_Surface_impl;

	struct Generate_Accessor : public Generate_Base {

		std::string absolutePhaseNamespace(Phase::T p) {
			return "SurfaceFrom" + representation.reprClassName + "_" + Phase::namespaceName(p);
		}

		virtual Generate_abstract_spec_Base*	 toGenerate_abstract_spec_Base    () { return nullptr; }
		virtual Generate_abstract_inco_Base*	 toGenerate_abstract_inco_Base    () { return nullptr; }

		virtual Generate_Surface_spec*			 toGenerate_SurfaceFromMRIS_spec  () { return nullptr; }
		virtual Generate_Surface_impl*			 toGenerate_SurfaceFromMRIS_impl  () { return nullptr; }
		virtual Generate_Surface_inco*			 toGenerate_SurfaceFromMRIS_inco  () { return nullptr; }

		virtual string classMemberPrefix(string const & classId)     = 0;
		virtual void beginClass (Phase::T p, string const & classId) = 0;
		virtual void endClass	(Phase::T p, string const & classId) = 0;

		virtual void generateAccessorClasses(Phase::T p)
		{
			string currClassId;
			bool writeCommented;
			ColumnedOutput::T* colsPtr = nullptr;
			vector<Prop*> deferredComGeneral;
			string classNameDotDot;
			bool writers;

			auto maybeEndC = [&](bool needsSep = true) {
				auto & cols = *colsPtr;
				if (toGenerate_abstract_spec_Base()) cols << endC;
				else if (needsSep) cols << " ";
			};

			auto maybeJust = [&](Justify j) {
				return toGenerate_abstract_spec_Base() ? j : Just_none;
			};

			walkClasses(representation,
				[&](Representation::PropHow const & propHow) {
					return propHow.prop->accessorClassId;
				},
				[&](string const & classId) {
					beginClass(p, classId);
					currClassId = classId;
					writeCommented = false;
					colsPtr = new ColumnedOutput::T;
					deferredComGeneral.clear();
					classNameDotDot = classMemberPrefix(classId);
					writers = false;
				},
				[&](Prop& prop, How* how, bool write) {
					typedef void os;
					typedef void indent;
					auto & cols = *colsPtr;

					if (currClassId == "Vertex" && !writers && write) {
						auto & cols = *colsPtr;
						if (toGenerate_abstract_spec_Base())
							cols << "inline ";
						else
							cols << ignoreWidth;
						cols << "void";
						maybeEndC();
						cols << classMemberPrefix(currClassId) << maybeJust(Just_left) << "which_coords";
						maybeEndC(false);
						cols << ignoreWidth
							<< "(int which, float *x, float *y, float *z) const";
						maybeEndC();
						generate_which_coords(p, currClassId, cols, false);
						cols << endR;
					}
					writers = write;

					if (prop.isGeneralComment()) {
						if (toGenerate_abstract_spec_Base()) {
							deferredComGeneral.push_back(&prop);
						}
						return;
					}

					if (!prop.type) return;
					if (!how) return;
					if (how->isImplDetail()) return;

					bool canRead = (prop.firstReadingPhase <= p);
					bool canWrite = (prop.firstWritingPhase <= p && p <= prop.lastWritingPhase);
					if (p == Phase::AllM) canWrite = (prop.firstWritingPhase != Phase::end);

					if ((!write && !canRead)
					||   (write && !canWrite)
					   ) {

						// Throw away general comments that precede variables that are thrown away.
						//
						deferredComGeneral.clear();

						// Skip this one
						return;
					}

					for (auto deferred : deferredComGeneral) {
						cols << ignoreWidth << "// " << deferred->comment << endR;
					}
					deferredComGeneral.clear();

					if (write && !writeCommented) {
						writeCommented = true; cols << endR << "" << endR;
					}

					if (toGenerate_abstract_spec_Base())
						cols << "inline ";
					else
						cols << ignoreWidth;

					string fname = how->getId(prop);;
					if (write) {
						fname = how->setId(prop);
						cols << "void";
					} else {
						cols << how->retType(prop);
					}
					maybeEndC();
					
					cols << classNameDotDot << maybeJust(Just_left) << fname;
					maybeEndC(false);

					cols << "(";
					maybeEndC(false);
					auto s = how->indexArg(prop);
					if (s.size()) {
						cols << s;
					}
					bool hasSep = false;
					if (write && s.size()) {
						cols << ",";
						hasSep = true;
					}
					maybeEndC(write && hasSep);
					if (write) cols << maybeJust(Just_right) << how->setToArg(prop);
					maybeEndC(false);
					cols << ")";
					if (!write) cols << " const";
					maybeEndC();
					
					generateBeforeComment(cols, prop, write);

					const char* comment = "  //";
					if (prop.repeatedSize) {
						cols << comment << " size() is " << prop.repeatedSize->id;
						comment = ".  ";
					}

					if (prop.comment.size()) cols << comment << "  " << prop.comment;
					cols << endR;
					
					generateMoreLines(cols, currClassId, prop, *how, write);
				},

				[&](string const & classId) {
					colsPtr->append(tos(), depth);
					endClass(p, classId);
					delete colsPtr; colsPtr = nullptr;
				});
		}

        virtual void generateBeforeComment(ColumnedOutput::T& cols, bool const write) {
			cols << ";";
        }
        
		virtual void generateBeforeComment(ColumnedOutput::T& cols, Prop const & d, bool const write) {
		    generateBeforeComment(cols, write);
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, string const & classId, Prop const & prop, How const & how, bool const write) {
		}

        virtual void generate_which_coords(Phase::T p, string const & classId, ColumnedOutput::T& cols, bool const write) {
			generateBeforeComment(cols, write);
        }

		void generateNamespace(string const & parent, string const & name, Representation & representation) {
			indent() << "namespace " << name << " {" << endl;
			depth++;
			indent() << "typedef " << representation.reprClassName << " Representation;" << endl;

			for (int doModifiers = 0; doModifiers < 2; doModifiers++) {
				for (auto p = Phase::T(); p < Phase::end; p = Phase::T(p + 1)) {
					auto pns = Phase::namespaceName(p);

					int isModifier = (pns.back() == 'M');
					if (!isModifier != !doModifiers) continue;

					ofstream *child_os = nullptr;
					if (toGenerate_abstract_spec_Base()) {
						auto child_fnm = parent + "_" + pns + ".h";
						indent() << "#include \"" + child_fnm + "\"" << endl;
						child_os = new ofstream(createFilesWithin + child_fnm);
						os_push(child_os);
					} else {
						tos() << endl << endl;
					}

					indent() << "namespace " << pns << " {" << endl;
						generateAccessorClasses(p);
					indent() << "} // namespace " << pns << endl;

					if (child_os) {
						os_pop();
						delete child_os;
					}
				}
			}

			if (toGenerate_abstract_inco_Base()) {
                                indent() << endl;
				indent() << "struct Repr_Elt { " << endl;
                depth++;
                indent() << "bool operator==(Repr_Elt const & rhs) const { return repr == rhs.repr && idx == rhs.idx; }" << endl;
                indent() << "bool operator!=(Repr_Elt const & rhs) const { return repr != rhs.repr || idx != rhs.idx; }" << endl;
                depth--;
                indent() << "protected: " << endl;
				depth++;
				indent() << "Representation* repr; size_t idx; " << endl;
				indent() << "Repr_Elt() : repr(nullptr), idx(0) {}" << endl;
				indent() << "Repr_Elt(Representation* repr, size_t idx) : repr(repr), idx(idx) {}" << endl;
				indent() << "Repr_Elt(Repr_Elt const & src) : repr(src.repr), idx(src.idx) {}" << endl;
				tos() << endl;
				for (auto p = Phase::T(); p < Phase::end; p = Phase::T(p + 1)) {
					auto isFriend = [&](string const & classId) {
						indent() << "friend struct " << name << "::" << Phase::namespaceName(p) << "::" << classId << ";" << endl;
					};
					isFriend("Face");
					isFriend("Vertex");
					isFriend("Surface");
				}
				depth--;
				indent() << "};" << endl;
			}

			depth--;
			indent() << "} // " << name << endl;
		}

		Generate_Accessor(ostream & os, Representation & representation) : Generate_Base(os, representation)
		{
		}
	};

	struct Generate_abstract_inco_Base : public Generate_Accessor {
        virtual Generate_abstract_inco_Base*	 toGenerate_abstract_inco_Base    () { return this; }

		Generate_abstract_inco_Base(
			ostream & os,
            string const & parent,
			Representation & representation) : Generate_Accessor(os, representation)
		{
		}

		virtual string classMemberPrefix(string const & classId) { return ""; }
		virtual void   beginClass(Phase::T p, string const & classId) {}
		virtual void   endClass  (Phase::T p, string const & classId) {}

		virtual void generateAccessorClasses(Phase::T p) {
			walkClasses(representation,
				[&](Representation::PropHow const & propHow) { return propHow.prop->accessorClassId; },
				[&](string const & classId) {
					depth++;
					indent() << "struct " << classId << ";" << endl;
					depth--; },
				[&](Prop& prop, How* how, bool write) {},
				[&](string const & classId) {});
		}

	};

	struct Generate_abstract_spec_Base : public Generate_Accessor {

		virtual Generate_abstract_spec_Base* toGenerate_abstract_spec_Base() { return this; }

		Generate_abstract_spec_Base(
            ostream & os,
            string const & parent,
			Representation & representation) : Generate_Accessor(os, representation)
		{
		}

		virtual string classMemberPrefix(string const & classId) {
			return "";
		}

		virtual void beginClass(Phase::T p, string const & classId) {
            auto const phaseNamespaceId = Phase::namespaceName(p);
            bool isSurface = (classId == "Surface");

			indent() << "struct " << classId << " : public Repr_Elt {" << endl;
			depth++;
            {
                if (classId != "Surface") { indent() << "typedef " << phaseNamespaceId << "::Surface Surface;" << endl; }
                if (classId != "Face")    { indent() << "typedef " << phaseNamespaceId << "::Face    Face;"    << endl; }
                if (classId != "Vertex")  { indent() << "typedef " << phaseNamespaceId << "::Vertex  Vertex;"  << endl; }
            }
            
			ColumnedOutput::T cols;
			cols << "inline " << classId << endC << "(" << endC << ""                                            << endC << ");" << endR;
			cols << "inline " << classId << endC << "(" << endC << "" << classId << " const & src"               << endC << ");" << endR;
			cols << "inline " << classId << endC << "(" << endC << "Representation* representation" << (isSurface?"":", size_t idx") << endC << ");" << endR;

			bool const isModifier = (phaseNamespaceId.back() == 'M');

			for (auto pLater = Phase::T(p + 1); pLater < Phase::end; pLater = Phase::T(pLater + 1)) {
				auto pnsLater = Phase::namespaceName(pLater);
				if (isModifier && pLater != Phase::AllM) continue;
				cols << "inline " << classId << endC
					<< "(" << endC << pnsLater << "::" << classId << " const & src" << endC << ");"
					<< endR;
			}
            
            if (classId == "Face") {
                cols << "int fno"     << endC << "() const { return idx; }" << endR;
            }
            if (classId == "Vertex") {
                cols << "int vno"     << endC << "() const { return idx; }" << endR;
            }
            if (classId == "Surface") {
                generateVariousSurfaceMethods(p, cols);
            }
                    
			cols.append(tos(), depth);
			tos() << endl;
		}

        void generateVariousSurfaceMethods(Phase::T p, ColumnedOutput::T & cols) {
            if (p == Phase::XYZPositionM || p == Phase::DistortM || p == Phase::DistortM) {
                cols << "void freeDistsButNotOrig() { MRISfreeDistsButNotOrig(repr); }" << endR;
            }
        }
        
		virtual void endClass(Phase::T p, string const & classId) {
			// indent() << absolutePns(p) << "_" << c.id << " // implementation details" << endl;
			depth--;
			indent() << "}; // " << classId << endl << endl;
		}

	};

	struct Generate_abstract_impl_Base : public Generate_Accessor {

		Generate_abstract_impl_Base(
            ostream & os, 
            string const & parent,
			Representation & representation) : Generate_Accessor(os, representation) {
		}

		virtual string classMemberPrefix(string const & classId) {
			return classId + "::";
		}

		virtual void beginClass(Phase::T p, string const & classId) {
			ColumnedOutput::T cols;
            
            bool isSurface = (classId == "Surface");

			cols << classMemberPrefix(classId) << classId << endC << "(" << endC << "" << endC << ") {}" << endR;
			cols << classMemberPrefix(classId) << classId << endC << "(" << endC << "Representation* representation" << (isSurface?"":", size_t idx") << endC
                << ") : Repr_Elt(representation," << (isSurface  ? "0" : "idx" ) << ") {}" << endR;
			cols << classMemberPrefix(classId) << classId << endC << "(" << endC << "" << classId << " const & src" << endC << ") : Repr_Elt(src) {}" << endR;

			bool const isModifier = (Phase::namespaceName(p).back() == 'M');

			for (auto pLater = Phase::T(p + 1); pLater < Phase::end; pLater = Phase::T(pLater + 1)) {
				if (isModifier && pLater != Phase::AllM) continue;
				cols << classMemberPrefix(classId) << classId << endC << "(" << endC << Phase::namespaceName(pLater) << "::" << classId << " const & src" << endC << ") : Repr_Elt(src) {}" << endR;
			}
			cols.append(tos(),depth);
			tos() << endl;
		}

		virtual void endClass(Phase::T p, string const & classId) {
			tos() << endl << endl;
		}

		virtual void generateBeforeComment(ColumnedOutput::T& cols, Prop const & d, bool const write) {
			cols << "{";
		}

	};


	struct Generate_Surface_inco : public Generate_abstract_inco_Base {
		virtual Generate_Surface_inco* toGenerate_Surface_inco() { return this; }
		Generate_Surface_inco(
			ostream & os,
			string const & parent,
			Representation & representation) : Generate_abstract_inco_Base(os, parent, representation)
		{
			generateNamespace(parent, "SurfaceFrom" + representation.reprClassName, representation);
		}
	};

	struct Generate_Surface_spec : public Generate_abstract_spec_Base {
		virtual Generate_Surface_spec* toGenerate_Surface_spec() { return this; }

		Generate_Surface_spec(
			ostream & os,
			string const & parent,
			Representation & representation) : Generate_abstract_spec_Base(os, parent, representation)
		{
			generateNamespace(parent, "SurfaceFrom" + representation.reprClassName, representation);
		}

	};

	struct Generate_Surface_impl : public Generate_abstract_impl_Base {
		virtual Generate_Surface_impl* toGenerate_Surface_impl() { return this; }

		Generate_Surface_impl(
			ostream & os,
			string const & parent,
			Representation & representation) : Generate_abstract_impl_Base(os, parent, representation)
		{
			generateNamespace(parent, "SurfaceFrom" + representation.reprClassName, representation);
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, string const & classId, Prop const & prop, How const & how, bool const write) {
			bool const hasIndexArg = how.indexArg(prop).size();
			cols << ignoreWidth << Just_left << "    ";

			string getSetExpr = how.getSetExpr(how.reprNoArrow(),prop,write);

			if (!write) cols << "return ";
			cols << getSetExpr << ";";
			cols << endR;
			cols << ignoreWidth << Just_left << "}" << endR;
		}

		virtual void generate_which_coords(Phase::T p, string const & classId, ColumnedOutput::T& cols, bool const write) {
			if (representation.reprClassName == "MRIS") 
				generate_which_coords_MRIS(p, classId, cols, write);
			else
				generate_which_coords_other(p, classId, cols, write);
		}
		
		void generate_which_coords_MRIS(Phase::T p, string const & classId, ColumnedOutput::T& cols, bool const write) {
			cols << ignoreWidth << Just_left << "{" << endR;
            cols << ignoreWidth << Just_left << "    " << endC;
			cols << "MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);" << endR;
            cols << ignoreWidth << Just_left << "}" << endR;
        }
		void generate_which_coords_other(Phase::T p, string const & classId, ColumnedOutput::T& cols, bool const write) {
			cols << ignoreWidth << Just_left << "{" << endR;
			cols << ignoreWidth << Just_left << "    " << endR;
			// cols << "repr->vertexCoord2XYZ(idx, which, x, y, z);" << endR;

			cols << "#define CASE(WHICH, FIELD) \\" << endR;
			cols << "  case WHICH: \\" << endR;
			cols << "    *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \\" << endR;
			cols << "    break;" << endR;
			cols << "" << endR;
			cols << "  switch (which) {" << endR;

			for (auto & propHow : representation.implements) {
				auto const & prop = *propHow.prop;
				auto const   how  =  propHow.how;
				if (!how) continue;
				if (prop.which == "") continue;
				bool canRead = (prop.firstReadingPhase <= p);
				if (!canRead) continue;
				assert(prop.type->id == "float");
				cols << ("    CASE(" + prop.which + "," + prop.id.substr(0, prop.id.size() - 1) + ")") << endR;
			}

			cols << "    default:" << endR;
			cols << "      *x = *y = *z = 0.0;" << endR;
			cols << "      ErrorExit(ERROR_UNSUPPORTED, \"which_coords: unsupported which %d\", which);" << endR;
			cols << "      break;" << endR;
			cols << "  }" << endR;
			cols << "" << endR;
			cols << "#undef CASE" << endR;

			cols << ignoreWidth << Just_left << "}" << endR;
		}
	};

}

// Section 6: main							Invoke some of the generators
//
namespace Representations {
	using namespace Abstract_Representation;
	static void build(std::vector<Representation*> & representations);
}

void generate(std::vector<Abstract_Representation::Representation*> & representations)
{
	using namespace std;

	for (auto representationP : representations) {
		auto & representation = *representationP;
		
		{	ofstream os(createFilesWithin + representation.rootHeaderFile);
			Representation_Generators::Generate_usingHow(os, representation);
		}

		string const root_fnm = "mrisurf_SurfaceFrom" + representation.reprClassName + "_generated";
		ofstream os(createFilesWithin + root_fnm + ".h");
		os << "#pragma once" << endl;

		auto fnm_inco = root_fnm + "_prefix";
		os << "#include \"./" << fnm_inco << ".h\"" << endl;
		{
			ofstream os_inco(createFilesWithin + fnm_inco + ".h");
			Accessor_Generators::Generate_Surface_inco(os_inco, fnm_inco, representation);
		}

		{
			Accessor_Generators::Generate_Surface_spec(os, root_fnm, representation);
		}

		auto fnm_impl = root_fnm + "_suffix";
		os << "#include \"./" << fnm_impl << ".h\"" << endl;
		{
			ofstream os_impl(createFilesWithin + fnm_impl + ".h");
			Accessor_Generators::Generate_Surface_impl(os_impl, fnm_impl, representation);
		}
	}

	auto cmd = "dir " + createFilesWithin;
	std::cout << cmd << endl;
	::system(cmd.c_str());
}


int main(int argc, const char* *argv)
{
	// Build the data structure that describes all about the surface, vertex, face, etc.
	std::vector<Abstract_Representation::Representation*> representations;
	Representations::build(representations);

	// Walk the data structure to output the sources needed for a variety of different purposes
	//
    createFilesWithin = (argc > 1) ? argv[1] : "./tmp";
	generate(representations);

	// Done
	return 0;
}

// Section 7: Representations						Define the properties, what each representation supports, what each phase can access
//
namespace Representations {

	// MRIS support
	//
	struct How_FACE : public HowDirect { 
		How_FACE() { this->m_secondaryName = "face_type_"; } 
		virtual string storeId(Prop const & prop) const {
			return "faces[idx]."+HowDirect::storeId(prop);
		}
	};
	struct How_FACE_v : public How_FACE {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return How_FACE::getSetExpr(reprNoArrow, prop, write);
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->faces[idx].v[i] = to.idx";
		}
	};
	struct How_VT : public HowDirect { 
		How_VT() { this->m_secondaryName = "VERTEX_TOPOLOGY"; } 
		virtual string storeId(Prop const & prop) const {
			return "vertices_topology[idx]." + HowDirect::storeId(prop);
		}
	};
	struct How_VT_indexed : public How_VT {
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
	};
	struct How_VT_f : public How_VT_indexed {
		virtual string retType(Prop const & prop) const { return "Face"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return How_VT_indexed::getSetExpr(reprNoArrow, prop, write);
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->vertices_topology[idx].f[i] = to.idx";
		}
	};
	struct How_VT_n : public How_VT_indexed {
		virtual string retType(Prop const & prop) const { return "size_t"; }
	};
	struct How_VT_e : public How_VT_indexed {
		virtual string retType(Prop const & prop) const { return "int"; }
	};
	struct How_VT_v : public How_VT_indexed {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return How_VT_indexed::getSetExpr(reprNoArrow, prop, write);
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->vertices_topology[idx].v[i] = to.idx";
		}
	};

	struct How_V : public HowDirect {
		How_V() { this->m_secondaryName = "vertex_type_"; } 
		virtual string storeId(Prop const & prop) const {
			return "vertices[idx]." + HowDirect::storeId(prop);
		}
	};
	struct How_V_indexed : public How_V {
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
	};
	struct How_V_indexed_float : public How_V_indexed {
		virtual string retType(Prop const & prop) const { return "float"; }
	};

	struct How_MRIS : public HowDirect { How_MRIS() { this->m_secondaryName = "MRIS"; } };
	struct How_MRIS_hidden : public How_MRIS { How_MRIS_hidden() { m_isImplDetail = true; } };
	struct How_MRIS_vertex : public How_MRIS {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return "Vertex(" + reprNoArrow  + ", " + reprNoArrow + "->" + prop.id + " - " + reprNoArrow + "->vertices)";
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->"+ prop.id +" = " + reprNoArrow + "->vertices + to.idx";
		}
	};
	struct How_MRIS_indexed : public How_MRIS {
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
	};
	struct How_MRIS_indexed_vertices : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return "Vertex(" + reprNoArrow + ", i)";
			return "TBD";
		}
	};
	struct How_MRIS_indexed_faces : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "Face"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return "Face(" + reprNoArrow + ", i)";
			return "TBD";
		}
	};
	struct How_MRIS_indexed_MRI_EDGE : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "MRI_EDGE"; }
	};
	struct How_MRIS_indexed_MRI_CORNER : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "MRI_CORNER"; }
	};
	struct How_MRIS_indexed_FaceNormCacheEntry : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "FaceNormCacheEntry"; }
	};
	struct How_MRIS_indexed_FaceNormDeferredEntry : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "FaceNormDeferredEntry"; }
	};
	struct How_MRIS_indexed_STRIP : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "STRIP"; }
	};
	struct How_MRIS_indexed_float : public How_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "float"; }
	};


	// MRISPV support
	//
	struct HowPV_Face : public HowDirect {
		HowPV_Face() { }
		virtual string memberType(Prop const & prop) const {	// the type of the field in the classId
			return storeType(prop) + "*";
		}
		virtual string memberId(Prop const & prop) const {		// the id of field in the classId
			return "f_" + prop.id;								// it might need to be indexed by [idx]
		}
		virtual string storeType(Prop const & prop) const {
			return prop.type->id;
		}
		virtual string storeId(Prop const & prop) const {
			return memberId(prop) + "[idx]";
		}
	};
	struct HowPV_F_indexed : public HowPV_Face {
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
	};
	struct HowPV_Face_v : public HowPV_Face {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return HowPV_Face::getSetExpr(reprNoArrow, prop, write);
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->f_v[idx][i] = to.idx";
		}
	};
	struct HowPV_V : public HowDirect {
		HowPV_V() { }
		virtual string memberType(Prop const & prop) const {	// the type of the field in the classId
			return prop.type->id + "*";
		}
		virtual string memberId(Prop const & prop) const {		// the id of field in the classId
			return "v_" + prop.id;								// it might need to be indexed by [idx]
		}
		virtual string storeType(Prop const & prop) const {
			return prop.type->id;
		}
		virtual string storeId(Prop const & prop) const {
			return memberId(prop) + "[idx]";
		}
	};
	struct HowPV_V_indexed : public HowPV_V {
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
	};
	struct HowPV_V_f : public HowPV_V_indexed {
		virtual string retType(Prop const & prop) const { return "Face"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write)        return "Face(" + reprNoArrow + ", " + reprNoArrow + "->v_f[idx][i])";
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->v_f[idx][i] = to.idx";
		}
	};
	struct HowPV_V_n : public HowPV_V_indexed {
		virtual string retType(Prop const & prop) const { return "size_t"; }
	};
	struct HowPV_V_e : public HowPV_V_indexed {
		virtual string retType(Prop const & prop) const { return "int"; }
	};
	struct HowPV_V_v : public HowPV_V_indexed {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return HowPV_V_indexed::getSetExpr(reprNoArrow, prop, write);
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->v_v[idx][i] = to.idx";
		}
	};

	struct HowPV_V_indexed_float : public HowPV_V_indexed {
		virtual string retType(Prop const & prop) const { return "float"; }
	};

	struct HowPV_MRIS : public HowDirect { HowPV_MRIS() { this->m_secondaryName = ""; } };
	struct HowPV_MRIS_hidden : public HowPV_MRIS { HowPV_MRIS_hidden() { m_isImplDetail = true; } };
	struct HowPV_MRIS_vertex : public HowPV_MRIS {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return "Vertex(" + reprNoArrow + ", " + reprNoArrow + "->" + prop.id + " - " + reprNoArrow + "->vertices)";
			return "cheapAssert(" + reprNoArrow + " == to.repr); " + reprNoArrow + "->" + prop.id + " = " + reprNoArrow + "->vertices + to.idx";
		}
	};
	struct HowPV_MRIS_indexed : public HowPV_MRIS {
		virtual string indexArg(Prop const & prop) const { return "size_t i"; }
	};
	struct HowPV_MRIS_indexed_vertices : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "Vertex"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return "Vertex(" + reprNoArrow + ", i)";
			return "TBD";
		}
	};
	struct HowPV_MRIS_indexed_faces : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "Face"; }
		virtual string getSetExpr(string const & reprNoArrow, Prop const & prop, bool write) const {
			if (!write) return "Face(" + reprNoArrow + ", i)";
			return "TBD";
		}
	};
	struct HowPV_MRIS_indexed_MRI_EDGE : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "MRI_EDGE"; }
	};
	struct HowPV_MRIS_indexed_MRI_CORNER : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "MRI_CORNER"; }
	};
	struct HowPV_MRIS_indexed_FaceNormCacheEntry : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "FaceNormCacheEntry"; }
	};
	struct HowPV_MRIS_indexed_FaceNormDeferredEntry : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "FaceNormDeferredEntry"; }
	};
	struct HowPV_MRIS_indexed_STRIP : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "STRIP"; }
	};
	struct HowPV_MRIS_indexed_float : public HowPV_MRIS_indexed {
		virtual string retType(Prop const & prop) const { return "float"; }
	};

	// MRIS_MP support
	//		It only supports a limited set of properties itself
	//		and an even smaller set are requested from an underlying MRIS or another MRIS_MP
	//
	//		It looks a lot like a subset of MRISPV
	//
	struct HowMP_likePV : public HowRedirect {
		HowMP_likePV(How* how) : HowRedirect(how) {}
		virtual string secondaryName() const { return ""; }
	};

	struct HowMP_fromMRIS : public HowRedirect {
		HowMP_fromMRIS(How * how) : HowRedirect(how) {}
		virtual string secondaryName() const { return ""; }
		virtual string reprNoArrow() const { return HowRedirect::reprNoArrow() + "->underlyingMRIS"; }
	};

	struct HowMP_fromMRIS_MP : public HowRedirect {
		HowMP_fromMRIS_MP(How* how) : HowRedirect(how) {}
		virtual string secondaryName() const { return ""; }
		virtual string reprNoArrow() const { return HowRedirect::reprNoArrow() + "->in_src"; }
	};

    struct HowMP_FaceNorm : public HowPV_Face {
        virtual string retType(Prop const & prop) const { return "FloatXYZ"; }
        virtual string storeType(Prop const & prop) const { return retType(prop); }
    };
	// Apply
	//
	struct RepresentationX : public Representation {
		RepresentationX(const char* rootHeaderFile, string const & surfaceClassName) : Representation(rootHeaderFile, surfaceClassName) { hows.push_back(nullptr); }
		std::vector<HowDirect*> hows;
	};

	static void doMRIS_MP(RepresentationX & rep_MRIS, RepresentationX & rep_MRISPV, RepresentationX & rep_MRIS_MP);

	static void build(std::vector<Representation*> & final_representations)
	{
		// Types are used to describe how properties are passed in (to Set), stored, and returned 
		//
		// Atomic types are passed in by value, stored by value, returned by value
		//
		auto t_bool						= new AtomicType("bool");

		auto t_char						= new AtomicType("char");
		auto t_short					= new AtomicType("short");
		auto t_int						= new AtomicType("int");
		auto t_long						= new AtomicType("long");
		auto t_uchar					= new AtomicType("uchar");
		auto t_ushort					= new AtomicType("ushort");
		auto t_uint						= new AtomicType("uint");
		auto t_ulong					= new AtomicType("ulong");

		auto t_double					= new AtomicType("double");
		auto t_float					= new AtomicType("float");
		auto t_constFloat				= new AtomicType("const float");

		auto t_pVoid					= new AtomicType("p_void");
		auto t_ppVoid					= new AtomicType("p_p_void");
		auto t_vertices_per_face_t		= new AtomicType("vertices_per_face_t");
		auto t_MATRIX					= new AtomicType("MATRIX");
		auto t_DMATRIX					= new AtomicType("DMATRIX");
		auto t_A3PDMATRIX				= new AtomicType("A3PDMATRIX");
		auto t_angles_per_triangle_t	= new AtomicType("angles_per_triangle_t");
		auto t_VOL_GEOM					= new AtomicType("VOL_GEOM");
		auto t_MRIS_cmdlines_t			= new AtomicType("MRIS_cmdlines_t");
		auto t_MRI						= new AtomicType("MRI");
		auto t_VERTEX					= new AtomicType("VERTEX");
		auto t_VERTEX_TOPOLOGY			= new AtomicType("VERTEX_TOPOLOGY");
		auto t_FACE						= new AtomicType("FACE");
		auto t_MRI_EDGE					= new AtomicType("MRI_EDGE");
		auto t_MRI_CORNER			    = new AtomicType("MRI_CORNER");
		auto t_FaceNormCacheEntry		= new AtomicType("FaceNormCacheEntry");
		auto t_FaceNormDeferredEntry	= new AtomicType("FaceNormDeferredEntry");
		auto t_STRIP					= new AtomicType("STRIP");
		auto t_LTA						= new AtomicType("LTA");
		auto t_MRIS_fname_t				= new AtomicType("MRIS_fname_t");
		auto t_MRIS_Status				= new AtomicType("MRIS_Status");
		auto t_MRIS_AREA_LABEL			= new AtomicType("MRIS_AREA_LABEL");
		auto t_MRIS_subject_name_t		= new AtomicType("MRIS_subject_name_t");
		auto t_COLOR_TABLE				= new AtomicType("COLOR_TABLE");

		// Pointer types are passed in by value, stored by value, returned by value
		//
		auto t_PMATRIX					= new PointerType("PMATRIX",			t_MATRIX);
		auto t_PDMATRIX					= new PointerType("PDMATRIX",			t_DMATRIX);
		auto t_PVERTEX					= new PointerType("PVERTEX",			t_VERTEX);
		auto t_PFACE					= new PointerType("PFACE",  			t_FACE);
		auto t_PLTA						= new PointerType("PLTA",				t_LTA);
		auto t_PMRIS_AREA_LABEL			= new PointerType("PMRIS_AREA_LABEL",	t_MRIS_AREA_LABEL);
		auto t_PCOLOR_TABLE				= new PointerType("PCOLOR_TABLE",		t_COLOR_TABLE);
		auto t_PMRI						= new PointerType("PMRI",				t_MRI);

		// These are passed in by TBD, stored as a count and a pointer, by TBD.
		//
		auto t_PR_float					= new PointerToRepeatedAtomicType("pSeveralFloat"					, t_float);
        auto t_PR_constFloat			= new PointerToRepeatedAtomicType("pSeveralConstFloat"              , t_constFloat);
		auto t_PR_int					= new PointerToRepeatedAtomicType("pSeveralInt"						, t_int);
		auto t_PR_uchar					= new PointerToRepeatedAtomicType("pSeveralUchar"					, t_uchar);
		auto t_PR_VERTEX				= new PointerToRepeatedAtomicType("pSeveralVERTEX"					, t_VERTEX);
		auto t_PR_VERTEX_TOPOLOGY		= new PointerToRepeatedAtomicType("pSeveralVERTEX_TOPOLOGY"			, t_VERTEX_TOPOLOGY);
		auto t_PR_FACE					= new PointerToRepeatedAtomicType("pSeveralFACE"					, t_FACE);
		auto t_PR_MRI_EDGE				= new PointerToRepeatedAtomicType("pSeveralMRI_EDGE"				, t_MRI_EDGE);
		auto t_PR_MRI_CORNER			= new PointerToRepeatedAtomicType("pSeveralMRI_CORNER"				, t_MRI_CORNER);
		auto t_PR_FaceNormCacheEntry	= new PointerToRepeatedAtomicType("pSeveralFaceNormCacheEntry"		, t_FaceNormCacheEntry);
		auto t_PR_FaceNormDeferredEntry = new PointerToRepeatedAtomicType("pSeveralFaceNormDeferredEntry"	, t_FaceNormDeferredEntry);
		auto t_PR_STRIP					= new PointerToRepeatedAtomicType("pSeveralSTRIP"					, t_STRIP);

		// As the following is executed, the properties go into all the representations currently in the set.
		//
		std::set<RepresentationX*> representations;

		string accessorClassId;

		auto rep_MRIS    = new RepresentationX("mrisurf_FACE_VERTEX_MRIS_generated.h"  , "MRIS");
		auto rep_MRISPV  = new RepresentationX("mrisurf_MRIS_PropertiesInVectors.h"    , "MRISPV");

		representations.insert(rep_MRIS);
		representations.insert(rep_MRISPV);

		auto howPush = [&](HowDirect* h0, HowDirect* h1) {
			rep_MRIS   ->hows.push_back(h0);
			rep_MRISPV ->hows.push_back(h1);
		};
		auto howMod = [&](HowDirect* h0, HowDirect* h1) {
			rep_MRIS   ->hows.back() = h0;
			rep_MRISPV ->hows.back() = h1;
		};
		auto howPop = [&]() {
			rep_MRIS   ->hows.pop_back();
			rep_MRISPV ->hows.pop_back();
		};

		auto makeProp = [&](Type* t, Id i, Phase::T rb, Phase::T wb, Phase::T we, string const & com = "") {
			if (rb == Phase::end && wb != Phase::end) rb = wb;
			auto added = new Prop(accessorClassId, t, i, rb, wb, we, com);
			for (auto & r : representations) {
				if (!r) continue;
				if (auto how = r->hows.back()) r->implements.push_back(Representation::PropHow(added, how));
			};
			return added;
		};

		auto addProp            = [&](Type* t, Id i, string const & com = "") { return makeProp(t, i, phaseRbegin, phaseWbegin, phaseWend, com); };
		auto addPropCom			= [&](string const & com) { return makeProp(nullptr, "", Phase::end, Phase::end, Phase::end, com); };
		auto addPropList		= [&](string const & com) { return addPropCom(com)->setCommentNature(ComList); };
		auto addPropListSublist = [&](string const & com) { return addPropCom(com)->setCommentNature(ComListSublist); };


		// Face
		accessorClassId = "Face";

		howPush(new How_FACE,new HowPV_Face);

		// Describe how the next few properties are to be stored
			phaseRbegin = phaseWbegin = phaseWend = Phase::TopologyM;
			addPropList("LIST_OF_FACE_ELTS");

			howPush(new How_FACE_v, new HowPV_Face_v);
			addProp(t_vertices_per_face_t, "v");
			howPop();

			phaseRbegin = phaseWbegin = phaseWend = Phase::XYZPositionConsequencesM;

		addProp(t_float, "area");
		addProp(t_angles_per_triangle_t, "angle");

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		addProp(t_angles_per_triangle_t, "orig_angle");

			phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::end;

		addProp(t_char, "ripflag");
		addProp(t_char, "oripflag");
		addProp(t_int,  "marked");

			phaseRbegin = phaseWbegin = phaseWend = Phase::XYZPositionConsequencesM;

		addProp(t_PDMATRIX, "norm");
		addProp(t_A3PDMATRIX, "gradNorm")->setNoHash();

		// Vertices
		accessorClassId = "Vertex";
		// First the VERTEX_TOPOLOGY members
			howPush(new How_VT,new HowPV_V);

			phaseRbegin = phaseWbegin = Phase::TopologyM; phaseWend = Phase::end;

			addPropList("LIST_OF_VERTEX_TOPOLOGY_ELTS");
			addPropCom("put the pointers before the ints, before the shorts, before uchars, to reduce size");
			addPropCom("the whole fits in much less than one cache line, so further ordering is no use");

			howPush(new How_VT_f, new HowPV_V_f);
	auto vtx_f =
		addProp(t_PR_int,		"f"                 , "array[v->num] the fno's of the neighboring faces         ");
			howMod(new How_VT_n, new HowPV_V_n);
	auto vtx_n =
		addProp(t_PR_uchar,		"n"            	    , "array[v->num] the face.v[*] index for this vertex        ");
			howMod(new How_VT_e, new HowPV_V_e);
		addProp(t_PR_int,		"e"                 , "edge state for neighboring vertices                      ");
			howPop();

			phaseRbegin = phaseWbegin = Phase::TopologyM; phaseWend = Phase::TopologyM;

			howPush(new How_VT_v,new HowPV_V_v);
	auto vtx_v =
		addProp(t_PR_int,		"v"                 , "array[v->vtotal or more] of vno, head sorted by hops     ");
			howPop();
	auto vtx_vnum = 								    
		addProp(t_short,		"vnum"              , "number of 1-hop neighbors    should use [p]VERTEXvnum(i) ");
													    
		addProp(t_short,		"v2num"             , "number of 1, or 2-hop neighbors                          ");
		addProp(t_short,		"v3num"             , "number of 1,2,or 3-hop neighbors                         ");
	auto vtx_vtotal =
		addProp(t_short,	    "vtotal"            , "total # of neighbors. copy of vnum.nsizeCur              ");
													    
		addProp(t_short,		"nsizeMaxClock"     , "copy of mris->nsizeMaxClock when v#num                   ")
			->setNoHash();							    
		addProp(t_uchar,		"nsizeMax"          , "the max nsize that was used to fill in vnum etc          ");
		addProp(t_uchar,		"nsizeCur"          , "index of the current v#num in vtotal                     ");
	auto vtx_num = 									    
		addProp(t_uchar,		"num"               , "number of neighboring faces                              ");

		vtx_f->setPRSize(vtx_num);
		vtx_n->setPRSize(vtx_num);
		vtx_v->setPRSize(vtx_vtotal);

		// The the VERTEX members
			howPush(new How_V, new HowPV_V);
	
			phaseRbegin = Phase::XYZPositionM; phaseWbegin = Phase::end; phaseWend = Phase::end;

			addPropList("LIST_OF_VERTEX_ELTS_1");
			addPropCom("managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]");

			howPush(new How_V_indexed_float, new HowPV_V_indexed_float);
	auto vtx_dist =
		addProp(t_PR_float,	"dist"			    , "distance to neighboring vertices based on  xyz   ");
			vtx_dist->setPRSize(vtx_vtotal);
	auto vtx_dist_orig =
		addProp(t_PR_float,	"dist_orig"		    , "distance to neighboring vertices based on origxyz");
			vtx_dist_orig->setPRSize(vtx_vtotal);
			howPop();

		addProp(t_int,      "dist_capacity"     , "-- should contain at least vtx_vtotal elements   ");
		addProp(t_int,      "dist_orig_capacity", "-- should contain at least vtx_vtotal elements   ");

			phaseRbegin = phaseWbegin = Phase::XYZPositionM; phaseWend = Phase::XYZPositionM;
			addPropCom("");

		addProp(t_float,	"x"					, "current coordinates	")->setWhich("CURRENT_VERTICES");
		addProp(t_float,	"y"					, "use MRISsetXYZ() to set");
		addProp(t_float,	"z"					);

			phaseRbegin = Phase::XYZPositionM; phaseWbegin = Phase::end; phaseWend = Phase::end;

		addPropCom("");
		addProp(t_float,	"origx"				, "original coordinates, see also MRIS::origxyz_status")->setWhich("ORIGINAL_VERTICES");
		addProp(t_float,	"origy"				, "use MRISsetOriginalXYZ(, ");
		addProp(t_float,	"origz"				, "or MRISsetOriginalXYZfromXYZ to set");

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		addPropCom("");
		addProp(t_float,	"nx")->setWhich("VERTEX_NORMALS");
		addProp(t_float,	"ny");
		addProp(t_float,	"nz", "curr normal");

			phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		addProp(t_float,	"pnx")->setWhich("PIAL_NORMALS");
		addProp(t_float,	"pny");
		addProp(t_float,	"pnz", "pial normal");

			addPropCom("");

		addProp(t_float,	"wnx")->setWhich("WHITE_NORMALS"); 
		addProp(t_float,	"wny"); 
		addProp(t_float,	"wnz", "white normal");
		addProp(t_float,	"onx"); 
		addProp(t_float,	"ony"); 
		addProp(t_float,	"onz", "original normal");
		addProp(t_float,	"dx"); 
		addProp(t_float,	"dy"); 
		addProp(t_float,	"dz", "current change in position");
		addProp(t_float,	"odx"); 
		addProp(t_float,	"ody"); 
		addProp(t_float,	"odz", "last change of position (for momentum, ");
		addProp(t_float,	"tdx"); 
		addProp(t_float,	"tdy"); 
		addProp(t_float,	"tdz", "temporary storage for averaging gradient");
		addProp(t_float,	"curv", "curr curvature");
		addProp(t_float,	"curvbak"); 
		addProp(t_float,	"val", "scalar data value (file: rh.val, sig2-rh.w)");
		addProp(t_float,	"imag_val", "imaginary part of complex data value");

			phaseRbegin = phaseWbegin = Phase::XYZPositionM; phaseWend = Phase::end;

		addProp(t_float,	"cx")->setWhich("CANONICAL_VERTICES");
		addProp(t_float,	"cy"); 
		addProp(t_float,	"cz", "coordinates in canonical coordinate system");

			phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		addProp(t_float,	"tx")->setWhich("TMP_VERTICES");
		addProp(t_float,	"ty"); 
		addProp(t_float,	"tz", "tmp coordinate storage");
		addProp(t_float,	"t2x")->setWhich("TMP2_VERTICES"); 
		addProp(t_float,	"t2y"); 
		addProp(t_float,	"t2z", "another tmp coordinate storage");
		addProp(t_float,	"targx"); 
		addProp(t_float,	"targy"); 
		addProp(t_float,	"targz", "target coordinates");
		addProp(t_float,	"pialx")->setWhich("PIAL_VERTICES"); 
		addProp(t_float,	"pialy"); 
		addProp(t_float,	"pialz", "pial surface coordinates");
		addProp(t_float,	"whitex")->setWhich("WHITE_VERTICES"); 
		addProp(t_float,	"whitey"); 
		addProp(t_float,	"whitez", "white surface coordinates");
		addProp(t_float,	"l4x"); 
		addProp(t_float,	"l4y"); 
		addProp(t_float,	"l4z", "layerIV surface coordinates");
		addProp(t_float,	"infx")->setWhich("INFLATED_VERTICES"); 
		addProp(t_float,	"infy"); 
		addProp(t_float,	"infz", "inflated coordinates");
		addProp(t_float,	"fx")->setWhich("FLATTENED_VERTICES"); 
		addProp(t_float,	"fy"); 
		addProp(t_float,	"fz", "flattened coordinates");
		addProp(t_int,		"px"); 
		addProp(t_int,		"qx"); 
		addProp(t_int,		"py"); 
		addProp(t_int,		"qy"); 
		addProp(t_int,		"pz"); 
		addProp(t_int,		"qz", "rational coordinates for exact calculations");
		addProp(t_float,	"e1x"); 
		addProp(t_float,	"e1y"); 
		addProp(t_float,	"e1z", "1st basis vector for the local tangent plane");
		addProp(t_float,	"e2x"); 
		addProp(t_float,	"e2y"); 
		addProp(t_float,	"e2z", "2nd basis vector for the local tangent plane");
		addProp(t_float,	"pe1x"); 
		addProp(t_float,	"pe1y"); 
		addProp(t_float,	"pe1z", "1st basis vector for the local tangent plane");
		addProp(t_float,	"pe2x"); 
		addProp(t_float,	"pe2y"); 
		addProp(t_float,	"pe2z", "2nd basis vector for the local tangent plane");

			addPropList("LIST_OF_VERTEX_ELTS_3");

        addProp(t_float,	"nc",		"curr length normal comp ");
        addProp(t_float,	"val2",		"complex comp data value (file: sig3-rh.w) ");
        addProp(t_float,	"valbak",	"scalar data stack ");
        addProp(t_float,	"val2bak",	"complex comp data stack ");
        addProp(t_float,	"stat",		"statistic ");
        addPropCom("");

        addProp(t_int,		"undefval",				"[previously dist=0] ");
        addProp(t_int,		"old_undefval",			"for smooth_val_sparse ");
        addProp(t_int,		"fixedval",				"[previously val=0] ");
        addPropCom("");

        addProp(t_float,	"fieldsign",			"fieldsign--final: -1, \"0\", \"1\" (file: rh.fs) ");
        addProp(t_float,	"fsmask",				"significance mask (file: rh.fm) ");
        addProp(t_float,	"d",					"for distance calculations ");

			addPropList("LIST_OF_VERTEX_ELTS_5");	

        addProp(t_int,		"annotation",			"area label (defunct--now from label file name!) ");
        addProp(t_char,		"oripflag");
        addProp(t_char,		"origripflag",			"cuts flags ");

			addPropList("LIST_OF_VERTEX_ELTS_7");

	    addProp(t_pVoid,	"vp",					"to store user's information ")->setNoHash();
        addProp(t_float,	"theta");
        addProp(t_float,	"phi",					"parameterization ");
		
			phaseRbegin = phaseWbegin = phaseWend = Phase::XYZPositionConsequencesM;
		
		addProp(t_float,	"area");

			phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		addProp(t_float,	"origarea");
        addProp(t_float,	"group_avg_area");
        addProp(t_float,	"K",					"Gaussian curvature ");
        addProp(t_float,	"H",					"mean curvature ");
        addProp(t_float,	"k1");
        addProp(t_float,	"k2",					"the principal curvatures ");
        addProp(t_float,	"mean");
        addProp(t_float,	"mean_imag",			"imaginary part of complex statistic ");
        addProp(t_float,	"std_error");
        addProp(t_uint,		"flags");
        addProp(t_int,		"fno",					"face that this vertex is in ");
        addProp(t_int,		"cropped");
        addProp(t_short,	"marked",				"for a variety of uses ");
        addProp(t_short,	"marked2");
        addProp(t_short,	"marked3");
        addProp(t_char,		"neg",					"1 if the normal vector is inverted ");
        addProp(t_char,		"border",				"flag ");

			phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::end;

		addProp(t_char,		"ripflag",				"vertex no longer exists - placed last to load the next vertex into cache");

		addPropList("LIST_OF_VERTEX_ELTS");
		addPropListSublist("LIST_OF_VERTEX_ELTS_1");
		addPropListSublist("LIST_OF_VERTEX_ELTS_3");
		addPropListSublist("LIST_OF_VERTEX_ELTS_5");
		addPropListSublist("LIST_OF_VERTEX_ELTS_7");

		// Edge

		// Surface
		accessorClassId = "Surface";
			howPush(new How_MRIS, new HowPV_MRIS);

			phaseWbegin = phaseWend = Phase::TopologyM;

			addPropList("LIST_OF_MRIS_ELTS_1");

			phaseRbegin = Phase::TopologyM; phaseWbegin = phaseWend = Phase::end;
    
		addPropCom("Fields being maintained by specialist functions");
		addProp(t_int,						"nverticesFrozen",	        "# of vertices on surface is frozen");
		addProp(t_int,						"nvertices",		        "# of vertices on surface, change by calling MRISreallocVerticesAndFaces et al");
		addProp(t_int,						"nfaces",		            "# of faces on surface, change by calling MRISreallocVerticesAndFaces et al");
		addProp(t_bool,						"faceAttachmentDeferred",	"defer connecting faces to vertices for performance reasons");
		addProp(t_int,						"nedges",		            "# of edges on surface");
        addProp(t_int,                      "ncorners",                 "# of triangle corners");
		addProp(t_int,						"nstrips");

			howPush(new How_MRIS_hidden, new HowPV_MRIS_hidden);
		addProp(t_PR_VERTEX_TOPOLOGY,		"vertices_topology");
			howMod(new How_MRIS_indexed_vertices, new HowPV_MRIS_indexed_vertices);
		addProp(t_PR_VERTEX,				"vertices");
			howPop();
		addProp(t_ppVoid,					"dist_storage",			    "the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored")->setNoHash();
		addProp(t_ppVoid,					"dist_orig_storage",		"the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored")->setNoHash();
		addProp(t_int,						"tempsAssigned",	        "State of various temp fields that can be borrowed if not already in use");

			phaseRbegin = Phase::TopologyM; phaseWbegin = phaseWend = Phase::end;

			howPush(new How_MRIS_indexed_faces, new HowPV_MRIS_indexed_faces);
		addProp(t_PR_FACE,					"faces");
			howMod(new How_MRIS_indexed_MRI_EDGE, new HowPV_MRIS_indexed_MRI_EDGE);
		addProp(t_PR_MRI_EDGE,				"edges");
			howMod(new How_MRIS_indexed_MRI_CORNER, new HowPV_MRIS_indexed_MRI_CORNER);
        addProp(t_PR_MRI_CORNER,            "corners");
			howMod(new How_MRIS_indexed_FaceNormCacheEntry, new HowPV_MRIS_indexed_FaceNormCacheEntry);
		addProp(t_PR_FaceNormCacheEntry,	"faceNormCacheEntries");
			howMod(new How_MRIS_indexed_FaceNormDeferredEntry, new HowPV_MRIS_indexed_FaceNormDeferredEntry);
		addProp(t_PR_FaceNormDeferredEntry,	"faceNormDeferredEntries");
			howPop();

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;
			howPush(new How_MRIS_indexed_STRIP, new HowPV_MRIS_indexed_STRIP);
		addProp(t_PR_STRIP,					"strips");
			howPop();

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		addProp(t_float,					"xctr");
		addProp(t_float,					"yctr");
		addProp(t_float,					"zctr");
		addProp(t_float,					"xlo");
		addProp(t_float,					"ylo");
		addProp(t_float,					"zlo");
		addProp(t_float,					"xhi");
		addProp(t_float,					"yhi");
		addProp(t_float,					"zhi");
		addProp(t_float,					"x0", "center of spherical expansion");
		addProp(t_float,					"y0");
		addProp(t_float,					"z0");

			addPropCom("v_temporal_pole, v_frontal_pole, and v_occipital_pole don't appear to be used, and are unusual being pointers to vertices");
			howPush(new How_MRIS_vertex, nullptr);
		addProp(t_PVERTEX,					"v_temporal_pole");			// WEIRD THAT THESE ARE POINTERS TO VERTICES
		addProp(t_PVERTEX,					"v_frontal_pole");
		addProp(t_PVERTEX,					"v_occipital_pole");
			howPop();
			addPropCom("");

		addProp(t_float,					"max_curv");
		addProp(t_float,					"min_curv");
		addProp(t_float,					"total_area");
		addProp(t_double,					"avg_vertex_area");

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::XYZPositionConsequencesM;

		addProp(t_double,					"avg_vertex_dist",			"set by MRIScomputeAvgInterVertexDist");
		addProp(t_double,					"std_vertex_dist");
		addProp(t_float,					"orig_area");
		addProp(t_float,					"neg_area");
		addProp(t_float,					"neg_orig_area",			"amount of original surface in folds");
		
			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		addProp(t_int,						"zeros");
		addProp(t_int,						"hemisphere",			"which hemisphere");
		
			phaseRbegin = Phase::ExistenceM; phaseWbegin = Phase::end; phaseWend = Phase::end;

		addProp(t_int,						"initialized");

			addPropList("LIST_OF_MRIS_ELTS_3");

		addProp(t_PLTA,						"lta");
		addProp(t_PMATRIX,					"SRASToTalSRAS_");
		addProp(t_PMATRIX,					"TalSRASToSRAS_");
		addProp(t_int,						"free_transform");
		addProp(t_double,					"radius",				"radius (if status==MRIS_SPHERE)");
		addProp(t_float,					"a");
		addProp(t_float,					"b");
		addProp(t_float,					"c",					"ellipsoid parameters");

			phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::ExistenceM;

		addProp(t_MRIS_fname_t,				"fname",				"file it was originally loaded from");

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		addProp(t_float,					"Hmin",					"min mean curvature");
		addProp(t_float,					"Hmax",					"max mean curvature");
		addProp(t_float,					"Kmin",					"min Gaussian curvature");
		addProp(t_float,					"Kmax",					"max Gaussian curvature");
		addProp(t_double,					"Ktotal",				"total Gaussian curvature");

			phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::ExistenceM;

		addProp(t_MRIS_Status,				"status",				"type of surface (e.g. sphere,"" plane)");
		addProp(t_MRIS_Status,				"origxyz_status",		"type of surface (e.g. sphere, plane) that this origxyz were obtained from");
		addProp(t_int,						"patch",				"if a patch of the surface");

			phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		addProp(t_int,						"nlabels");
		addProp(t_PMRIS_AREA_LABEL,			"labels",				"nlabels of these (may be null)");
		
			phaseRbegin = phaseWbegin = Phase::end; phaseWend = Phase::end;

		addProp(t_char,						"nsize",				"size of neighborhoods or -1");
		addProp(t_uchar,					"vtotalsMightBeTooBig", "MRISsampleDistances sets this");
		addProp(t_short,					"nsizeMaxClock",		"changed whenever an edge is added or removed, which invalidates the vertex v#num values")->setNoHash();
		addProp(t_char,						"max_nsize",			"max the neighborhood size has been set to (typically 3)");
		addProp(t_char,						"dist_nsize",			"max mrisComputeVertexDistances has computed distances out to");
		addProp(t_char,						"dist_orig_nsize",		"max mrisComputeOriginalVertexDistances has computed distances out to");
		addProp(t_char,						"dist_alloced_flags",	"two flags, set when any dist(1) or dist_orig(2) allocated");
		addProp(t_float,					"avg_nbrs",				"mean # of vertex neighbors");

			phaseRbegin = phaseWbegin = Phase::XYZPositionM; phaseWend = Phase::end;

		addProp(t_pVoid,					"vp",					"for misc. use")->setNoHash();
		addProp(t_float,					"alpha",				"rotation around z-axis");
		addProp(t_float,					"beta",					"rotation around y-axis");
		addProp(t_float,					"gamma",				"rotation around x-axis");
		addProp(t_float,					"da");
		addProp(t_float,					"db");
		addProp(t_float,					"dg",					"old deltas");
		addProp(t_int,						"type",					"what type of surface was this initially");

			phaseRbegin = Phase::ExistenceM; phaseWbegin = Phase::end; phaseWend = Phase::end;

		addProp(t_int,						"max_vertices",					"may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces");
		addProp(t_int,						"max_faces",					"may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces");

		addProp(t_MRIS_subject_name_t,		"subject_name",					"name of the subject");
		addProp(t_float,					"canon_area");
		addProp(t_int,						"noscale",						"don't scale by surface area if true");
			howPush(new How_MRIS_indexed_float, new HowPV_MRIS_indexed_float);
		addProp(t_PR_float,					"dx2",							"an extra set of gradient (not always alloced)");
		addProp(t_PR_float,					"dy2");
		addProp(t_PR_float,					"dz2");
			howPop();
		addProp(t_PCOLOR_TABLE,				"ct");
		addProp(t_int,						"orig_xyzspace",       "xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space");
		addProp(t_int,						"useRealRAS",					"if 0 (default), vertex position is a conformed volume RAS with c_(r,\"a\",\"s\")=0.  "
																				"else is a real RAS (volume stored RAS)");
		addProp(t_VOL_GEOM,					"vg",							"volume info from which this surface is created. valid iff vg.valid = 1");
		addProp(t_MRIS_cmdlines_t,			"cmdlines")->setNoHash();
		addProp(t_int,						"ncmds");
		addProp(t_float,					"group_avg_surface_area",		"average of total surface area for group");
		addProp(t_int,						"group_avg_vtxarea_loaded",		"average vertex area for group at each vertex");
		addProp(t_int,						"triangle_links_removed",		"for quad surfaces");
		addProp(t_pVoid,					"user_parms",					"for whatever the user wants to hang here")->setNoHash();
		addProp(t_PMATRIX,					"m_sras2vox",					"for converting surface ras to voxel");
		addProp(t_PMRI,						"mri_sras2vox",					"volume that the above matrix is for");
		addProp(t_pVoid,					"mht")->setNoHash();
		addProp(t_pVoid,					"temps")->setNoHash();

			addPropList("LIST_OF_MRIS_ELTS");
			addPropListSublist("LIST_OF_MRIS_ELTS_1");
			addPropListSublist("LIST_OF_MRIS_ELTS_3");
		

		// rather than polluting the above list with the MRIS_MP information, that is pulled into here...
		//
		auto rep_MRIS_MP = new RepresentationX("mrisurf_MRIS_MPPropertiesInVectors.h", "MRIS_MP");
		doMRIS_MP(*rep_MRIS, *rep_MRISPV, *rep_MRIS_MP);

		final_representations.clear();
		final_representations.push_back(rep_MRIS);
		final_representations.push_back(rep_MRISPV);
		final_representations.push_back(rep_MRIS_MP);
	}

	static void doMRIS_MP(RepresentationX & rep_MRIS, RepresentationX & rep_MRISPV, RepresentationX & rep_MRIS_MP)
	{
        
		struct PropHowMap : public std::map<string, Representation::PropHow*> {
			PropHowMap(RepresentationX & r) {
				for (auto & ph : r.implements) (*this)[ph.prop->key()] = &ph;
			}
		} rep_MRISPV_propHowMap(rep_MRISPV);

		enum Action { NotImpl, Stored, FromMRIS };
		std::map<string,Action> actions;
		auto insert = [&](const char* id, Action action = Stored) {
			actions[id] = action;
		};
		{ 
			insert("Surface.status",				FromMRIS);
			insert("Surface.patch",				        FromMRIS);
			insert("Surface.noscale",			        FromMRIS);
			insert("Surface.origxyz_status");
			insert("Surface.nvertices");
			insert("Surface.vertices");
			insert("Surface.nfaces");
                        insert("Surface.faces");
                        insert("Surface.faceNormCacheEntries");
                        insert("Surface.faceNormDeferredEntries");
			insert("Surface.nsize");
			insert("Surface.radius");
			insert("Surface.faces_topology");
			insert("Surface.dist_nsize");
			insert("Surface.xctr");
			insert("Surface.yctr");
			insert("Surface.zctr");
			insert("Surface.xlo");
			insert("Surface.xhi");
			insert("Surface.ylo");
			insert("Surface.yhi");
			insert("Surface.zlo");
			insert("Surface.zhi");
			insert("Surface.total_area");
			insert("Surface.avg_vertex_area");
			insert("Surface.avg_vertex_dist");
			insert("Surface.std_vertex_dist");
			insert("Surface.orig_area");
			insert("Surface.neg_orig_area");
			insert("Surface.neg_area");

			insert("Vertex.ripflag");
			insert("Vertex.VSize");
			insert("Vertex.whitex");
			insert("Vertex.whitey");
			insert("Vertex.whitez");
			insert("Vertex.pialx");
			insert("Vertex.pialy");
			insert("Vertex.pialz");
			insert("Vertex.wnx");
			insert("Vertex.wny");
			insert("Vertex.wnz");
			insert("Vertex.dist_orig");
			insert("Vertex.x");
			insert("Vertex.y");
			insert("Vertex.z");
			insert("Vertex.dx");
			insert("Vertex.dy");
			insert("Vertex.dz");
			insert("Vertex.dist_capacity");
			insert("Vertex.border");
			insert("Vertex.cx");
			insert("Vertex.cy");
			insert("Vertex.cz");
			insert("Vertex.curv");
			insert("Vertex.origx");
			insert("Vertex.origy");
			insert("Vertex.origz");
			insert("Vertex.origarea");
			insert("Vertex.fno");
			insert("Vertex.area");
			insert("Vertex.nx");
			insert("Vertex.ny");
			insert("Vertex.nz");
			insert("Vertex.neg");
			insert("Vertex.dist");
			
			insert("Face.ripflag");
			insert("Face.v");
			insert("Face.norm_orig_area");
			insert("Face.orig_angle");
			insert("Face.area");
			insert("Face.normSet");
			insert("Face.norm");
			insert("Face.angle");
		}		

		auto const how_implementationDetail      = new HowDirect();
		how_implementationDetail->m_isImplDetail = true;
        
        auto implementationDetail = [&](
            string type, string id, string comment)
        {
		    rep_MRIS_MP.implements.push_back(
			    Representation::PropHow(
				    new Prop("",new AtomicType(type),id, Phase::end, Phase::ExistenceM,comment),
				    how_implementationDetail));
        };
        
		implementationDetail("MRIS*",                   "underlyingMRIS", "for properties that are read from the underlying MRIS");
		implementationDetail("MRIS_MP*",                "in_src",         "since the in are not written, they can be shared by copies");
		implementationDetail("int",                     "in_ref_count",   "check the src doesn't go away");
        implementationDetail("VERTEX_TOPOLOGY const *", "vertices_topology", "pointer copied from MRIS");
        implementationDetail("FACE_TOPOLOGY   const *", "faces_topology",    "pointer copied from MRIS");
		
		implementationDetail("int*",                    "v_VSize",          "");
		implementationDetail("float**",                 "v_dist_buffer",    "");
		implementationDetail("const float*",            "f_norm_orig_area", "");
		implementationDetail("char*",                   "f_normSet",        "");

		for (size_t i = 0; i < rep_MRIS.implements.size(); i++) {
			auto & propHow = rep_MRIS.implements[i];
			auto   prop = propHow.prop;
            auto   iPropKey = prop->key();
			auto it = actions.find(iPropKey);
			if (it == actions.end()) continue;
            
            
			switch (it->second) {
				case NotImpl:
					break;
				case FromMRIS:
					rep_MRIS_MP.implements.push_back(Representation::PropHow(prop, new HowMP_fromMRIS(propHow.how)));
					break;
				case Stored: {
					auto it = rep_MRISPV_propHowMap.find(iPropKey);
					assert(it != rep_MRISPV_propHowMap.end());
					auto & pvPropHowImpl = *it->second;
                    How* how = new HowMP_likePV(pvPropHowImpl.how);
                    if (iPropKey == "Face.norm") {
                        std::cout << "Found " << iPropKey << std::endl;
                        how = new HowMP_FaceNorm;
                    }
                    rep_MRIS_MP.implements.push_back(Representation::PropHow(prop, how));
				}	break;
			}
		}
	}

}
