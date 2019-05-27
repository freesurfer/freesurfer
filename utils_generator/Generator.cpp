// Compile with     g++ -std=c++11 Generator.cpp
// Run with         rm -rf tmp ; mkdir tmp ; ./a.out ./tmp/
// Merge with       bcompare ./tmp ../include
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

static std::ostream& indent(std::ostream & os, size_t depth) {
	for (size_t i = 0; i < depth; i++) os << "    ";
	return os;
}

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
				for (size_t ci = 0; ci < r.size(); ci++) {
					os << sep;
					if (ci < colSizes.size() && colSizes[ci] > 0) 
						os << setw(colSizes[ci]);
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


namespace Generator {
	using namespace std;
	using namespace ColumnedOutput;

	static string uppercase(string const & s_init) {
		string s = s_init;
		transform(s.begin(), s.end(), s.begin(), ::toupper);
		return s;
	}

	string createFilesWithin;

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

	struct Clss; Clss* newClss(Id i);
	struct Dmbr;
	struct Fmbr;
	struct Type;

	enum CommentNature { ComGeneral, ComList, ComListSublist };

	const char* currentMrisVector = nullptr;

	struct DmbrModifiable {
		const char *  mrisVector;
		CommentNature commentNature;
		bool  nohash;
		Dmbr* repeatedSize;	// type must be a PointerToRepeatedType type, this is the element that gives the number of repeats
		DmbrModifiable() : mrisVector(currentMrisVector), commentNature(ComGeneral), nohash(false), repeatedSize(nullptr) {}
	};
	struct Dmbr : public DmbrModifiable {
		Clss*		const clss;
		Type*		const type;
		Id			const id;
		string      const comment;
        Phase::T    const firstReadingPhase;
		Phase::T	const firstWritingPhase;
		Phase::T	const lastWritingPhase;
		Dmbr(Clss* c, Type* t, Id i,              Phase::T wb             , string const & com = "") : clss(c), type(t), id(i), firstReadingPhase(wb), firstWritingPhase(wb), lastWritingPhase(wb), comment(com) {}
		Dmbr(Clss* c, Type* t, Id i,              Phase::T wb, Phase::T we, string const & com = "") : clss(c), type(t), id(i), firstReadingPhase(wb), firstWritingPhase(wb), lastWritingPhase(we), comment(com) {}
		Dmbr(Clss* c, Type* t, Id i, Phase::T rb, Phase::T wb, Phase::T we, string const & com = "") : clss(c), type(t), id(i), firstReadingPhase(rb), firstWritingPhase(wb), lastWritingPhase(we), comment(com) {}

		bool isGeneralComment() const { return !type && commentNature==ComGeneral; }

		Dmbr* setCommentNature(CommentNature to) { commentNature = to;     return this; }
		Dmbr* setNoHash()		  { nohash = true;     return this; }
		Dmbr* setPRSize(Dmbr* to) { repeatedSize = to; return this; }
	};

	struct Fmbr {
		Clss*		const clss;
		Type*		const retType;
		Id			const id;
		Clss*		const args;
		bool		const isWriter;
		Fmbr(Clss* c, Type* rt, Id i, Phase::T w, bool isWriter = false) : clss(c), retType(retType), id(i), args(newClss("<args>")), isWriter(isWriter) {}
	};

	Phase::T phaseRbegin = Phase::end;
	Phase::T phaseWbegin = Phase::end;
	Phase::T phaseWend   = Phase::end;

	struct Clss {
		Clss(Id i) : id(i) {}

		Id const id;
		map<Id,Dmbr*> dmbrMap;
		vector<Dmbr*> dmbr;
		vector<Fmbr*> fmbr;

		Dmbr* addDmbrCom(string const & com)				  { return addDmbr(nullptr, "", Phase::end, Phase::end, Phase::end, com); }
		Dmbr* addDmbrList(string const & com)				  { return addDmbrCom(com)->setCommentNature(ComList); }
		Dmbr* addDmbrListSublist(string const & com)		  { return addDmbrCom(com)->setCommentNature(ComListSublist); }
		Dmbr* addDmbr(Type* t, Id i, string const & com = "") { return addDmbr(t, i, phaseRbegin, phaseWbegin, phaseWend, com); }
		Dmbr* addDmbr(Type* t, Id i, Phase::T rb, Phase::T wb, Phase::T we, string const & com = "") {
            if (rb == Phase::end && wb != Phase::end) rb = wb;
			auto res = new Dmbr(this, t, i, rb, wb, we, com); 
			dmbr.push_back(res);
			dmbrMap.insert(make_pair(id,res));
			return res; 
		}
		
		Dmbr* findDmbr(Id i) {
			auto it = dmbrMap.find(i);
			return (it == dmbrMap.end()) ? nullptr : it->second;
		}

	}; Clss* newClss(Id i) { return new Clss(i); }

	struct AtomicType; struct PointerType; struct ArrayOfPointerType; struct PointerToRepeatedType;
	struct Type {
		Id			const id;
		Type(Id i) : id(i) {}
		virtual AtomicType*  toAtomicType () { return nullptr; }
		virtual PointerType* toPointerType() { return nullptr; }
		virtual ArrayOfPointerType*  toArrayOfPointerType() { return nullptr; }
		virtual PointerToRepeatedType* toPointerToRepeatedType() { return nullptr; }
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
	
	struct PointerToRepeatedType : public PointerType {
		virtual PointerToRepeatedType* toPointerToRepeatedType() { return this; }
		PointerToRepeatedType(Id i, AtomicType* at) : PointerType(i, at) {}
	};

	typedef vector<Clss*> Built;
	static Built build();

	struct Generate_Base {
		ostream & tos() { return *os.back(); }
		Built const & built;
		size_t depth;
		ostream& indent() {
			return ::indent(tos(),depth);
		}
		Generate_Base(
			ostream & init_os, 
			Built const & built) : os(), built(built), depth(0) {
			os.push_back(&init_os);
			tos() << endl;
			indent() << "// GENERATED SOURCE - DO NOT DIRECTLY EDIT" << endl;
			indent() << "// " << endl;
			indent() << "// =======================================" << endl;
                        indent() << "#define SEPARATE_VERTEX_TOPOLOGY" << endl;

		}
		void os_push(ostream * new_tos) { os.push_back(new_tos); }
		void os_pop() { os.pop_back(); }
	private:
		vector<ostream *> os;
	};

	struct Generate_mris : public Generate_Base {
		void generateClass(Clss & c) {
			bool const isVertex = (c.id == "Vertex");

			for (int pass = 0; pass < (isVertex ? 2 : 1); pass++) {
				bool isInVertexTopology  = false;
				bool doingVertexTopology = (pass == 0);

				string traditional_className = c.id;
				if (isVertex && pass == 0) traditional_className = "VERTEX_TOPOLOGY";
				if (traditional_className == "Surface") traditional_className = "MRIS";
				if (traditional_className == "Face")    traditional_className = "face_type_";
				if (traditional_className == "Vertex")  traditional_className = "vertex_type_";

				indent() << "struct " 
					<< traditional_className
					<< " {" << endl;
				depth++;

				ColumnedOutput::T cols;

				for (auto & d : c.dmbr) {

					if (d->commentNature == ComList) 
						isInVertexTopology = d->comment == "LIST_OF_VERTEX_TOPOLOGY_ELTS";
					if (isVertex && (doingVertexTopology != isInVertexTopology)) continue;

					if (!d->type && d->commentNature != ComGeneral) continue;

					const char * comment = "//";
					if (d->type) {
						cols << d->type->id << endC << Just_left << d->id << endC << ";";
						comment = "  //";
						if (d->repeatedSize) {
							cols << ignoreWidth << comment << " size() is " << d->repeatedSize->id;
							comment = ".  ";
						}
					}
					if (d->comment.size()) cols << ignoreWidth << comment << "  " << d->comment;
					cols << endR;
				}
				cols.append(tos(), depth);

				depth--;
				indent() << "};		// " << traditional_className << endl << endl;
			}
		}

		void generateMacros(Clss & c) {

			bool eltEmitted  = false;
			bool endlPending = false;
			bool endMPending = false;
			auto emitEndOfMacroIfNeeded = [&] {
				if (!endMPending) return;
				depth--;
				eltEmitted = false;
				if (endlPending) tos() << " \\" << endl; endlPending = false;
				if (endMPending) indent() << "// end of macro" << endl << endl; endMPending = false;
			};

			bool isVertexTopology = false;

			for (auto & d : c.dmbr) {

				if (d->isGeneralComment()) continue;

				if (eltEmitted ) { tos() << " SEP";		eltEmitted  = false; }
				if (endlPending) { tos() << " \\" << endl; endlPending = false; }

				if (!d->type) {
					switch (d->commentNature) {
					case ComList:
						emitEndOfMacroIfNeeded();
						indent() << "#define " << d->comment << " \\" << endl; eltEmitted = false; endlPending = false; endMPending = true;
						depth++;
						break;
					case ComListSublist:
						indent() << d->comment;  eltEmitted = true; endlPending = true;
						break;
					}
				} else if (endMPending) {
					const char* macro = (d->nohash) ? "ELTX" : "ELTT";
					auto t = d->type;
					if (auto tp = t->toPointerType()) { macro = "ELTP"; t = tp->target; }
					indent() << macro << "(" << t->id << "," << d->id << ") "; eltEmitted = true; endlPending = true;
				}
			}

			emitEndOfMacroIfNeeded();
		}

		Generate_mris(
			ostream & os, 
			Built const & built) : Generate_Base(os, built)
		{
			for (auto & cp : built) generateClass(*cp);
			for (auto & cp : built) generateMacros(*cp);
		}
	};

	struct Generate_SurfaceFromMRIS_inco;
	struct Generate_SurfaceFromMRIS_spec;
	struct Generate_SurfaceFromMRIS_impl;

	struct Generate_abstract_Base : public Generate_Base {

		std::string absolutePns(Phase::T p) {
			return "SurfaceFromMRIS_" + Phase::namespaceName(p);
		}

		virtual Generate_SurfaceFromMRIS_inco* toGenerate_SurfaceFromMRIS_inco() { return nullptr; }
		virtual Generate_SurfaceFromMRIS_spec* toGenerate_SurfaceFromMRIS_spec() { return nullptr; }
		virtual Generate_SurfaceFromMRIS_impl* toGenerate_SurfaceFromMRIS_impl() { return nullptr; }

		virtual string classPrefix(Clss & c) = 0;
		virtual void   beginClass (Phase::T p, Clss & c) = 0;
		virtual void   endClass	  (Phase::T p, Clss & c) = 0;

		string baseClass(Phase::T p, Clss& c) {
#if 1
			return "MRIS_Elt";
#else
			string result;
			auto const enhances = Phase::enhances[p];
			if (enhances == Phase::end) {
				result = "MRIS_Elt";
			}
			else {
				result = Phase::namespaceName(enhances);
				result += "::";
				result += c.id;
			}
			return result;
#endif
		}

		Id retTypeId(Clss const & c, Dmbr const & d) {
			auto t = d.type;
			auto tr = t->toPointerToRepeatedType();

			if (c.id == "Vertex") {
				if (d.id == "f") return "Face";
				if (d.id == "n") return "size_t";
				if (d.id == "v") return "Vertex";

			} else if (c.id == "Face") {
				if (d.id == "v") return "Vertex";

			} else if (c.id == "Surface") {
			    if (d.id == "vertices") return "Vertex";
			    if (d.id == "faces")    return "Face";
			}

			if (tr) return tr->target->id;

            if (t->id == "PVERTEX") return "Vertex";
            if (t->id == "PFACE")   return "Face";
            
			return t->id;
		}

		virtual void generateClass(Phase::T p, Clss & c) {
			beginClass(p, c);

			for (int w = 0; w < 2; w++) {
				auto const write = (w==1);
				
				bool writeCommented = false;

				ColumnedOutput::T cols; {
					typedef void os;
					typedef void indent;

					auto maybeEndC = [&](bool needsSep = true) {
						if (toGenerate_SurfaceFromMRIS_spec()) cols << endC; 
						else if (needsSep) cols << " ";
					};

					auto maybeJust = [&](Justify j) {
						return toGenerate_SurfaceFromMRIS_spec() ? j : Just_none;
					};
					
					vector<Dmbr*> deferredComGeneral;

					for (auto & dp : c.dmbr) {
						auto const & d = *dp;

						if (d.isGeneralComment()) {
							if (toGenerate_SurfaceFromMRIS_spec()) {
								deferredComGeneral.push_back(dp);
							}
							continue;
						}

						if (!d.type) continue;
					
                        bool canRead  = (d.firstReadingPhase <= p);
                        bool canWrite = (d.firstWritingPhase <= p && p <= d.lastWritingPhase);
						if (p == Phase::AllM) canWrite = (d.firstWritingPhase != Phase::end);
                            
						if ((!write && !canRead)
						||  ( write && !canWrite)
                           ) {

							// Throw away general comments that precede variables that are thrown away.
							//
							deferredComGeneral.clear();

                            // Skip this one
							continue;
						}

                        if (c.id == "Surface") {
                            if (d.id == "vertices_topology") continue;
                        }
                        
						for (auto deferred : deferredComGeneral) {
							cols << ignoreWidth << "// " << deferred->comment << endR;
						}
						deferredComGeneral.clear();

						if (write && !writeCommented) { 
							writeCommented = true; cols << endR << "" << endR;
						}

						auto t   = d.type;
						auto tr  = t->toPointerToRepeatedType();
						auto rti = retTypeId(c, d);

						if (toGenerate_SurfaceFromMRIS_spec()) 
							cols << "inline ";
						else
							cols << ignoreWidth;

						string fname = d.id;
						if (write) {
							cols << "void";
							fname = "set_" + d.id;
						} else {
							cols << rti;
						}
						maybeEndC();
						cols << classPrefix(c) << maybeJust(Just_left) << fname;
						maybeEndC(false);
						
						cols << "(";
						maybeEndC(false);
						bool needsSep = false;
                        bool indexingAVector = (c.id=="Face" && d.id=="v");
						if (t->toPointerToRepeatedType() || indexingAVector) {
							cols << "size_t i";
							needsSep = true;
						} 
						if (write && (tr || indexingAVector)) {
							cols << ",";
							needsSep = true;
						}
						maybeEndC(write && needsSep);
						if (write) cols << maybeJust(Just_right) << rti << " to";
						maybeEndC(false);
						cols << ")";
						if (!write) cols << " const";
						maybeEndC();
						generateBeforeComment(cols,d,write);
						
						const char* comment = "  //";
						if (tr && d.repeatedSize) {
							cols << comment << " size() is " << d.repeatedSize->id;
							comment = ".  ";
						}
						if (d.comment.size()) cols << comment << "  " << d.comment;
						cols << endR;
						generateMoreLines(cols,c,d,write);
					}

				} cols.append(tos(),depth);
			}
			endClass(p, c);
		}

		virtual void generateBeforeComment(ColumnedOutput::T& cols, Dmbr const & d, bool const write) {
			cols << ";";
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, Clss const & c, Dmbr const & d, bool const write) {
		}

		virtual void generateNamespaceMacros(Phase::T p) {
		}

		void generateNamespace(string const & parent) {
			indent() << "namespace SurfaceFromMRIS {" << endl;
			depth++;

			for (int doModifiers = 0; doModifiers < 2; doModifiers++) {
				for (auto p = Phase::T(); p < Phase::end; p = Phase::T(p + 1)) {
					auto pns = Phase::namespaceName(p);

					int isModifier = (pns.back() == 'M');
					if (!isModifier != !doModifiers) continue;

					ofstream *child_os = nullptr;
					if (toGenerate_SurfaceFromMRIS_spec()) {
						auto child_fnm = parent + "_" + pns + ".h";
						indent() << "#include \"" + child_fnm + "\"" << endl;
						child_os = new ofstream(createFilesWithin + child_fnm);
						os_push(child_os);
					} else {
						tos() << endl << endl;
					}

					indent() << "namespace " << pns << " {" << endl;
						generateNamespaceMacros(p);
						for (auto & cp : built) generateClass(p, *cp);
					indent() << "} // namespace " << pns << endl;

					if (child_os) {
						os_pop();
						delete child_os;
					}
				}
			}

			if (toGenerate_SurfaceFromMRIS_inco()) {
				indent() << "struct MRIS_Elt { " << endl;
                depth++;
                indent() << "bool operator==(MRIS_Elt const & rhs) const { return mris == rhs.mris && idx == rhs.idx; }" << endl;
                indent() << "bool operator!=(MRIS_Elt const & rhs) const { return mris != rhs.mris || idx != rhs.idx; }" << endl;
                depth--;
                indent() << "protected: " << endl;
				depth++;
				indent() << "MRIS* mris; size_t idx; " << endl;
				indent() << "MRIS_Elt() : mris(nullptr), idx(0) {}" << endl;
				indent() << "MRIS_Elt(MRIS* mris, size_t idx) : mris(mris), idx(idx) {}" << endl;
				indent() << "MRIS_Elt(MRIS_Elt const & src) : mris(src.mris), idx(src.idx) {}" << endl;
				tos() << endl;
				for (auto p = Phase::T(); p < Phase::end; p = Phase::T(p + 1)) {
					for (auto & cp : built) {
						indent() << "friend struct SurfaceFromMRIS::" << Phase::namespaceName(p) << "::" << cp->id << ";" << endl;
					}
				}
				depth--;
				indent() << "};" << endl;
			}

			depth--;
			indent() << "} // namespace SurfaceFromMRIS" << endl;
		}

		Generate_abstract_Base(ostream & os, Built const & built) : Generate_Base(os, built)
		{
		}
	};

	struct Generate_SurfaceFromMRIS_inco : public Generate_abstract_Base {
		virtual Generate_SurfaceFromMRIS_inco* toGenerate_SurfaceFromMRIS_inco() { return this; }

		Generate_SurfaceFromMRIS_inco(
			ostream & os,
            string const & parent,
			Built const & built) : Generate_abstract_Base(os, built)
		{
			generateNamespace(parent);
		}

		virtual void generateNamespaceMacros(Phase::T p) {
		}

		virtual string classPrefix(Clss & c) { return ""; }
		virtual void   beginClass(Phase::T p, Clss & c) {}
		virtual void   endClass(Phase::T p, Clss & c) {}

		virtual void generateClass(Phase::T p, Clss & c) {
			depth++;
			indent() << "struct " << c.id << ";" << endl;
			depth--;
		}
	};

	struct Generate_SurfaceFromMRIS_spec : public Generate_abstract_Base {
		virtual Generate_SurfaceFromMRIS_spec* toGenerate_SurfaceFromMRIS_spec() { return this; }

		Generate_SurfaceFromMRIS_spec(
            ostream & os,
            string const & parent,
            Built const & built) : Generate_abstract_Base(os, built)
		{
			generateNamespace(parent);
		}

		virtual string classPrefix(Clss & c) {
			return "";
		}

		virtual void beginClass(Phase::T p, Clss & c) {
            bool isSurface = (c.id == "Surface");

			indent() << "struct " << c.id << " : public " << baseClass(p,c) << " {" << endl;
			depth++;

			ColumnedOutput::T cols;
			cols << "inline " << c.id << endC << "(" << endC << ""                                            << endC << ");" << endR;
			cols << "inline " << c.id << endC << "(" << endC << "" << c.id << " const & src"                  << endC << ");" << endR;
			cols << "inline " << c.id << endC << "(" << endC << "MRIS* mris" << (isSurface?"":", size_t idx") << endC << ");" << endR;

			bool const isModifier = (Phase::namespaceName(p).back() == 'M');

			for (auto pLater = Phase::T(p + 1); pLater < Phase::end; pLater = Phase::T(pLater + 1)) {
				auto pnsLater = Phase::namespaceName(pLater);
				if (isModifier && pLater != Phase::AllM) continue;
				cols << "inline " << c.id << endC
					<< "(" << endC << pnsLater << "::" << c.id << " const & src" << endC << ");" 
					<< endR;
			}
            
            if (c.id == "Face") {
                cols << "int fno"     << endC << "() const { return idx; }" << endR;
            }
            if (c.id == "Vertex") {
                cols << "int vno"     << endC << "() const { return idx; }" << endR;
            }
            
			cols.append(tos(), depth);
			tos() << endl;
		}

		virtual void endClass(Phase::T p, Clss & c) {
			// indent() << absolutePns(p) << "_" << c.id << " // implementation details" << endl;
			depth--;
			indent() << "};" << endl << endl;
		}

	};

	struct Generate_SurfaceFromMRIS_impl : public Generate_abstract_Base {
		virtual Generate_SurfaceFromMRIS_impl* toGenerate_SurfaceFromMRIS_impl() { return this; }

		Generate_SurfaceFromMRIS_impl(
            ostream & os, 
            string const & parent,
            Built const & built) : Generate_abstract_Base(os, built) {
			generateNamespace(parent);
		}

		virtual string classPrefix(Clss & c) {
			return c.id + "::";
		}

		virtual void beginClass(Phase::T p, Clss & c) {
			ColumnedOutput::T cols;
            
            bool isSurface = (c.id == "Surface");

			cols << classPrefix(c) << c.id << endC << "(" << endC << "" << endC << ") {}" << endR;
			cols << classPrefix(c) << c.id << endC << "(" << endC << "MRIS* mris" << (isSurface?"":", size_t idx") << endC 
                << ") : " << baseClass(p, c) << "(mris," << (isSurface  ? "0" : "idx" ) << ") {}" << endR;
			cols << classPrefix(c) << c.id << endC << "(" << endC << "" << c.id << " const & src" << endC << ") : " << baseClass(p, c) << "(src) {}" << endR;

			bool const isModifier = (Phase::namespaceName(p).back() == 'M');

			for (auto pLater = Phase::T(p + 1); pLater < Phase::end; pLater = Phase::T(pLater + 1)) {
				if (isModifier && pLater != Phase::AllM) continue;
				cols << classPrefix(c) << c.id << endC << "(" << endC << Phase::namespaceName(pLater) << "::" << c.id << " const & src" << endC << ") : " << baseClass(p, c) << "(src) {}" << endR;
			}
			cols.append(tos(),depth);
			tos() << endl;
		}

		virtual void endClass(Phase::T p, Clss & c) {
			tos() << endl << endl;
		}

		virtual void generateBeforeComment(ColumnedOutput::T& cols, Dmbr const & d, bool const write) {
			cols << "{";
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, Clss const & c, Dmbr const & d, bool const write) {
			auto t = d.type;
			auto tr = t->toPointerToRepeatedType();
			cols << ignoreWidth << Just_left << "    ";

				string raw = "mris->";
				if (d.mrisVector) raw += d.mrisVector, raw += "[idx].";
				raw += d.id;
				if (tr) raw += "[i]";

				if (c.id == "Face" && d.id == "v") {
					if (!write) {
						cols << "return Vertex(mris,mris->faces[idx].v[i])";
					} else {
						cols << "cheapAssert(mris == to.mris); mris->faces[idx].v[i] = to.idx";
					}
				} else if (c.id == "Vertex" && d.id == "f") {
					if (!write) {
						cols << "return Face(mris," << raw << ")";
					} else {
						cols << "cheapAssert(mris == to.mris); " << raw << " = to.idx";
					}
				} else if (c.id == "Vertex" && d.id == "v") {
					if (!write) {
						cols << "return Vertex(mris," << raw << ")";
					}
					else {
						cols << "cheapAssert(mris == to.mris); " << raw << " = to.idx";
					}
				} else if (c.id == "Surface" && (d.id == "v_temporal_pole" || d.id == "v_frontal_pole" || d.id == "v_occipital_pole" )) {
					if (!write) {
						cols << "return Vertex(mris," << raw << " - mris->vertices)";
					}
					else {
						cols << "cheapAssert(mris == to.mris); " << raw << " = mris->vertices + to.idx";
					}
				} else if (c.id == "Surface" && d.id == "vertices" ) {
					if (!write) {
						cols << "return Vertex(mris,i)";
					}
					else {
						assert(false);
					}
				} else if (c.id == "Surface" && d.id == "faces" ) {
					if (!write) {
						cols << "return Face(mris,i)";
					}
					else {
						assert(false);
					}
				} else {
					string cvtPrefix = "", cvtFromPrefix = "";
					string cvtSuffix = "", cvtFromSuffix = "";
					if (c.id == "Vertex" && d.id == "n") {
						cvtPrefix = "size_t("; cvtFromPrefix = tr->target->id + "(";
						cvtSuffix = cvtFromSuffix = ")";
					}

					if (!write) cols << "return " << cvtPrefix;
					cols << raw;
					if (write) cols << " = " << cvtFromPrefix << "to" << cvtFromSuffix; else cols << cvtSuffix;
				}
				cols << ";" << endR;
			cols << ignoreWidth << Just_left << "}" << endR;
		}
	};

	void generate(Built const & built) {
		if (true) {
			ofstream os(createFilesWithin + "mrisurf_FACE_VERTEX_MRI_traditional.h");
			Generate_mris(os, built);
		}
		if (true) {
			string const root_fnm = "mrisurf_SurfaceFromMRIS_generated";
			ofstream os(createFilesWithin + root_fnm + ".h");

			auto fnm_inco = root_fnm + "_prefix";
			os << "#include \"./" << fnm_inco << ".h\"" << endl;
			{
				ofstream os_inco(createFilesWithin + fnm_inco + ".h");
				Generate_SurfaceFromMRIS_inco(os_inco, fnm_inco, built);
			}

			{
				Generate_SurfaceFromMRIS_spec(os, root_fnm, built);
			}

			auto fnm_impl = root_fnm + "_suffix";
			os << "#include \"./" << fnm_impl << ".h\"" << endl;
			{
				ofstream os_impl(createFilesWithin + fnm_impl + ".h");
				Generate_SurfaceFromMRIS_impl(os_impl, fnm_impl, built);
			}
		}
		system(("dir " + createFilesWithin).c_str());
	}
}


int main(int argc, const char* *argv)
{
	// Build the data structure that describes all about the surface, vertex, face, etc.
	auto clssVec = Generator::build();

	// Walk the data structure to output the sources needed for a variety of different purposes
	//
    Generator::createFilesWithin = (argc > 1) ? argv[1] : "./tmp";
	Generator::generate(clssVec);

	// Done
	return 0;
}


namespace Generator {

	static Built build() {

		auto t_double	= new AtomicType("double");
		auto t_float	= new AtomicType("float");
		auto t_char		= new AtomicType("char");
		auto t_short	= new AtomicType("short");
		auto t_int		= new AtomicType("int");
		auto t_long		= new AtomicType("long");
		auto t_uchar	= new AtomicType("uchar");
		auto t_ushort	= new AtomicType("ushort");
		auto t_uint		= new AtomicType("uint");
		auto t_ulong	= new AtomicType("ulong");
		auto t_bool		= new AtomicType("bool");

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
		auto t_FaceNormCacheEntry		= new AtomicType("FaceNormCacheEntry");
		auto t_FaceNormDeferredEntry	= new AtomicType("FaceNormDeferredEntry");
		auto t_STRIP					= new AtomicType("STRIP");
		auto t_LTA						= new AtomicType("LTA");
		auto t_MRIS_fname_t				= new AtomicType("MRIS_fname_t");
		auto t_MRIS_Status				= new AtomicType("MRIS_Status");
		auto t_MRIS_AREA_LABEL			= new AtomicType("MRIS_AREA_LABEL");
		auto t_MRIS_subject_name_t		= new AtomicType("MRIS_subject_name_t");
		auto t_COLOR_TABLE				= new AtomicType("COLOR_TABLE");

		auto t_PMATRIX					= new PointerType("PMATRIX",			t_MATRIX);
		auto t_PDMATRIX					= new PointerType("PDMATRIX",			t_DMATRIX);
		auto t_PVERTEX					= new PointerType("PVERTEX",			t_VERTEX);
		auto t_PFACE					= new PointerType("PFACE",  			t_FACE);
		auto t_PLTA						= new PointerType("PLTA",				t_LTA);
		auto t_PMRIS_AREA_LABEL			= new PointerType("PMRIS_AREA_LABEL",	t_MRIS_AREA_LABEL);
		auto t_PCOLOR_TABLE				= new PointerType("PCOLOR_TABLE",		t_COLOR_TABLE);
		auto t_PMRI						= new PointerType("PMRI",				t_MRI);

		auto t_PR_float					= new PointerToRepeatedType("pSeveralFloat"					, t_float);
		auto t_PR_int					= new PointerToRepeatedType("pSeveralInt"					, t_int);
		auto t_PR_uchar					= new PointerToRepeatedType("pSeveralUchar"					, t_uchar);
		auto t_PR_VERTEX				= new PointerToRepeatedType("pSeveralVERTEX"				, t_VERTEX);
		auto t_PR_VERTEX_TOPOLOGY		= new PointerToRepeatedType("pSeveralVERTEX_TOPOLOGY"		, t_VERTEX_TOPOLOGY);
		auto t_PR_FACE					= new PointerToRepeatedType("pSeveralFACE"					, t_FACE);
		auto t_PR_MRI_EDGE				= new PointerToRepeatedType("pSeveralMRI_EDGE"				, t_MRI_EDGE);
		auto t_PR_FaceNormCacheEntry	= new PointerToRepeatedType("pSeveralFaceNormCacheEntry"	, t_FaceNormCacheEntry);
		auto t_PR_FaceNormDeferredEntry = new PointerToRepeatedType("pSeveralFaceNormDeferredEntry"	, t_FaceNormDeferredEntry);
		auto t_PR_STRIP					= new PointerToRepeatedType("pSeveralSTRIP"					, t_STRIP);

		auto vtx = new Clss("Vertex");
		auto fac = new Clss("Face");
		// auto edg = new Clss("Edge");
		auto sur = new Clss("Surface");

		currentMrisVector = "faces";
		phaseRbegin = phaseWbegin = phaseWend = Phase::TopologyM;

		fac->addDmbrList("LIST_OF_FACE_ELTS");

		fac->addDmbr(t_vertices_per_face_t, "v");

		phaseRbegin = phaseWbegin = phaseWend = Phase::XYZPositionConsequencesM;

		fac->addDmbr(t_float, "area");
		fac->addDmbr(t_angles_per_triangle_t, "angle");

		phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		fac->addDmbr(t_angles_per_triangle_t, "orig_angle");

		phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::end;

		fac->addDmbr(t_char, "ripflag");
		fac->addDmbr(t_char, "oripflag");
		fac->addDmbr(t_int, "marked");

		phaseRbegin = phaseWbegin = phaseWend = Phase::XYZPositionConsequencesM;

		fac->addDmbr(t_PDMATRIX, "norm");
		fac->addDmbr(t_A3PDMATRIX, "gradNorm")->setNoHash();

		currentMrisVector = "vertices_topology";
		phaseRbegin = phaseWbegin = Phase::TopologyM; phaseWend = Phase::end;

		// vtx->dmbr();
		//	LIST_OF_VERTEX_TOPOLOGY_ELTS
		vtx->addDmbrList("LIST_OF_VERTEX_TOPOLOGY_ELTS");
		vtx->addDmbrCom("put the pointers before the ints, before the shorts, before uchars, to reduce size");
		vtx->addDmbrCom("the whole fits in much less than one cache line, so further ordering is no use");

	auto vtx_f =
		vtx->addDmbr(t_PR_int,		"f"                 , "array[v->num] the fno's of the neighboring faces         ");
	auto vtx_n =
		vtx->addDmbr(t_PR_uchar,	"n"            	    , "array[v->num] the face.v[*] index for this vertex        ");
		vtx->addDmbr(t_PR_int,		"e"                 , "edge state for neighboring vertices                      ");

		phaseRbegin = phaseWbegin = Phase::TopologyM; phaseWend = Phase::TopologyM;

	auto vtx_v =
		vtx->addDmbr(t_PR_int,		"v"                 , "array[v->vtotal or more] of vno, head sorted by hops     ");
	auto vtx_vnum = 								    
		vtx->addDmbr(t_short,		"vnum"              , "number of 1-hop neighbors    should use [p]VERTEXvnum(i, ");
													    
		vtx->addDmbr(t_short,		"v2num"             , "number of 1, or 2-hop neighbors                          ");
		vtx->addDmbr(t_short,		"v3num"             , "number of 1,2,or 3-hop neighbors                         ");
	auto vtx_vtotal =
		vtx->addDmbr(t_short,	    "vtotal"            , "total # of neighbors. copy of vnum.nsizeCur              ");
													    
		vtx->addDmbr(t_short,		"nsizeMaxClock"     , "copy of mris->nsizeMaxClock when v#num                   ")
			->setNoHash();							    
		vtx->addDmbr(t_uchar,		"nsizeMax"          , "the max nsize that was used to fill in vnum etc          ");
		vtx->addDmbr(t_uchar,		"nsizeCur"          , "index of the current v#num in vtotal                     ");
	auto vtx_num = 									    
		vtx->addDmbr(t_uchar,		"num"               , "number of neighboring faces                              ");

		vtx_f->setPRSize(vtx_num);
		vtx_n->setPRSize(vtx_num);
		vtx_v->setPRSize(vtx_vtotal);

		currentMrisVector = "vertices";
	
    phaseRbegin = Phase::XYZPositionM; phaseWbegin = Phase::end; phaseWend = Phase::end;

		vtx->addDmbrList("LIST_OF_VERTEX_ELTS_1");

		vtx->addDmbrCom("managed by MRISfreeDists[_orig] and MRISmakeDists[_orig]");
	auto vtx_dist =
		vtx->addDmbr(t_PR_float,	"dist"			    , "distance to neighboring vertices based on  xyz   ");
	auto vtx_dist_orig =									
		vtx->addDmbr(t_PR_float,	"dist_orig"		    , "distance to neighboring vertices based on origxyz");
		vtx->addDmbr(t_int,         "dist_capacity"     , "-- should contain at least vtx_vtotal elements   ");
		vtx->addDmbr(t_int,         "dist_orig_capacity", "-- should contain at least vtx_vtotal elements   ");

		vtx_dist->setPRSize     (vtx_vtotal);
		vtx_dist_orig->setPRSize(vtx_vtotal);

		phaseRbegin = phaseWbegin = Phase::XYZPositionM; phaseWend = Phase::XYZPositionM;

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"x"					, "current coordinates	");
		vtx->addDmbr(t_float,		"y"					, "use MRISsetXYZ() to set");
		vtx->addDmbr(t_float,		"z"					);

		phaseRbegin = Phase::XYZPositionM; phaseWbegin = Phase::end; phaseWend = Phase::end;

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"origx"				, "original coordinates, see also MRIS::origxyz_status");
		vtx->addDmbr(t_float,		"origy"				, "use MRISsetOriginalXYZ(, ");
		vtx->addDmbr(t_float,		"origz"				, "or MRISsetOriginalXYZfromXYZ to set");

		phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"nx");
		vtx->addDmbr(t_float,		"ny");
		vtx->addDmbr(t_float,		"nz", "curr normal");

		phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		vtx->addDmbr(t_float,		"pnx");
		vtx->addDmbr(t_float,		"pny");
		vtx->addDmbr(t_float,		"pnz", "pial normal");

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"wnx"); 
		vtx->addDmbr(t_float,		"wny"); 
		vtx->addDmbr(t_float,		"wnz", "white normal");
		vtx->addDmbr(t_float,		"onx"); 
		vtx->addDmbr(t_float,		"ony"); 
		vtx->addDmbr(t_float,		"onz", "original normal");
		vtx->addDmbr(t_float,		"dx"); 
		vtx->addDmbr(t_float,		"dy"); 
		vtx->addDmbr(t_float,		"dz", "current change in position");
		vtx->addDmbr(t_float,		"odx"); 
		vtx->addDmbr(t_float,		"ody"); 
		vtx->addDmbr(t_float,		"odz", "last change of position (for momentum, ");
		vtx->addDmbr(t_float,		"tdx"); 
		vtx->addDmbr(t_float,		"tdy"); 
		vtx->addDmbr(t_float,		"tdz", "temporary storage for averaging gradient");
		vtx->addDmbr(t_float,		"curv", "curr curvature");
		vtx->addDmbr(t_float,		"curvbak"); 
		vtx->addDmbr(t_float,		"val", "scalar data value (file: rh.val, sig2-rh.w)");
		vtx->addDmbr(t_float,		"imag_val", "imaginary part of complex data value");

		phaseRbegin = phaseWbegin = Phase::XYZPositionM; phaseWend = Phase::end;

		vtx->addDmbr(t_float,		"cx");
		vtx->addDmbr(t_float,		"cy"); 
		vtx->addDmbr(t_float,		"cz", "coordinates in canonical coordinate system");

		phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		vtx->addDmbr(t_float,		"tx");
		vtx->addDmbr(t_float,		"ty"); 
		vtx->addDmbr(t_float,		"tz", "tmp coordinate storage");
		vtx->addDmbr(t_float,		"tx2"); 
		vtx->addDmbr(t_float,		"ty2"); 
		vtx->addDmbr(t_float,		"tz2", "tmp coordinate storage");
		vtx->addDmbr(t_float,		"targx"); 
		vtx->addDmbr(t_float,		"targy"); 
		vtx->addDmbr(t_float,		"targz", "target coordinates");
		vtx->addDmbr(t_float,		"pialx"); 
		vtx->addDmbr(t_float,		"pialy"); 
		vtx->addDmbr(t_float,		"pialz", "pial surface coordinates");
		vtx->addDmbr(t_float,		"whitex"); 
		vtx->addDmbr(t_float,		"whitey"); 
		vtx->addDmbr(t_float,		"whitez", "white surface coordinates");
		vtx->addDmbr(t_float,		"l4x"); 
		vtx->addDmbr(t_float,		"l4y"); 
		vtx->addDmbr(t_float,		"l4z", "layerIV surface coordinates");
		vtx->addDmbr(t_float,		"infx"); 
		vtx->addDmbr(t_float,		"infy"); 
		vtx->addDmbr(t_float,		"infz", "inflated coordinates");
		vtx->addDmbr(t_float,		"fx"); 
		vtx->addDmbr(t_float,		"fy"); 
		vtx->addDmbr(t_float,		"fz", "flattened coordinates");
		vtx->addDmbr(t_int,			"px"); 
		vtx->addDmbr(t_int,			"qx"); 
		vtx->addDmbr(t_int,			"py"); 
		vtx->addDmbr(t_int,			"qy"); 
		vtx->addDmbr(t_int,			"pz"); 
		vtx->addDmbr(t_int,			"qz", "rational coordinates for exact calculations");
		vtx->addDmbr(t_float,		"e1x"); 
		vtx->addDmbr(t_float,		"e1y"); 
		vtx->addDmbr(t_float,		"e1z", "1st basis vector for the local tangent plane");
		vtx->addDmbr(t_float,		"e2x"); 
		vtx->addDmbr(t_float,		"e2y"); 
		vtx->addDmbr(t_float,		"e2z", "2nd basis vector for the local tangent plane");
		vtx->addDmbr(t_float,		"pe1x"); 
		vtx->addDmbr(t_float,		"pe1y"); 
		vtx->addDmbr(t_float,		"pe1z", "1st basis vector for the local tangent plane");
		vtx->addDmbr(t_float,		"pe2x"); 
		vtx->addDmbr(t_float,		"pe2y"); 
		vtx->addDmbr(t_float,		"pe2z", "2nd basis vector for the local tangent plane");

		vtx->addDmbrList("LIST_OF_VERTEX_ELTS_3");

        vtx->addDmbr(t_float,		"nc",		"curr length normal comp ");
        vtx->addDmbr(t_float,		"val2",		"complex comp data value (file: sig3-rh.w) ");
        vtx->addDmbr(t_float,		"valbak",	"scalar data stack ");
        vtx->addDmbr(t_float,		"val2bak",	"complex comp data stack ");
        vtx->addDmbr(t_float,		"stat",		"statistic ");
        vtx->addDmbrCom("");

        vtx->addDmbr(t_int,			"undefval",				"[previously dist=0] ");
        vtx->addDmbr(t_int,			"old_undefval",			"for smooth_val_sparse ");
        vtx->addDmbr(t_int,			"fixedval",				"[previously val=0] ");
        vtx->addDmbrCom("");

        vtx->addDmbr(t_float,		"fieldsign",			"fieldsign--final: -1, \"0\", \"1\" (file: rh.fs) ");
        vtx->addDmbr(t_float,		"fsmask",				"significance mask (file: rh.fm) ");
        vtx->addDmbr(t_float,		"d",					"for distance calculations ");

	    vtx->addDmbrList("LIST_OF_VERTEX_ELTS_5");	

        vtx->addDmbr(t_int,			"annotation",			"area label (defunct--now from label file name!) ");
        vtx->addDmbr(t_char,		"oripflag");
        vtx->addDmbr(t_char,		"origripflag",			"cuts flags ");

	    vtx->addDmbrList("LIST_OF_VERTEX_ELTS_7");

	    vtx->addDmbr(t_pVoid,		"vp",					"to store user's information ")->setNoHash();
        vtx->addDmbr(t_float,		"theta");
        vtx->addDmbr(t_float,		"phi",					"parameterization ");
		
	phaseRbegin = phaseWbegin = phaseWend = Phase::XYZPositionConsequencesM;
		
		vtx->addDmbr(t_float,		"area");

	phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		vtx->addDmbr(t_float,		"origarea");
        vtx->addDmbr(t_float,		"group_avg_area");
        vtx->addDmbr(t_float,		"K",					"Gaussian curvature ");
        vtx->addDmbr(t_float,		"H",					"mean curvature ");
        vtx->addDmbr(t_float,		"k1");
        vtx->addDmbr(t_float,		"k2",					"the principal curvatures ");
        vtx->addDmbr(t_float,		"mean");
        vtx->addDmbr(t_float,		"mean_imag",			"imaginary part of complex statistic ");
        vtx->addDmbr(t_float,		"std_error");
        vtx->addDmbr(t_uint,		"flags");
        vtx->addDmbr(t_int,			"fno",					"face that this vertex is in ");
        vtx->addDmbr(t_int,			"cropped");
        vtx->addDmbr(t_short,		"marked",				"for a variety of uses ");
        vtx->addDmbr(t_short,		"marked2");
        vtx->addDmbr(t_short,		"marked3");
        vtx->addDmbr(t_char,		"neg",					"1 if the normal vector is inverted ");
        vtx->addDmbr(t_char,		"border",				"flag ");

	phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::end;

		vtx->addDmbr(t_char,		"ripflag",				"vertex no longer exists - placed last to load the next vertex into cache");

		vtx->addDmbrList("LIST_OF_VERTEX_ELTS");
		vtx->addDmbrListSublist("LIST_OF_VERTEX_ELTS_1");
		vtx->addDmbrListSublist("LIST_OF_VERTEX_ELTS_3");
		vtx->addDmbrListSublist("LIST_OF_VERTEX_ELTS_5");
		vtx->addDmbrListSublist("LIST_OF_VERTEX_ELTS_7");

		// edg->dmbr();

		currentMrisVector = nullptr;
		phaseWbegin = phaseWend = Phase::TopologyM;

		sur->addDmbrList("LIST_OF_MRIS_ELTS_1");

	phaseRbegin = Phase::TopologyM; phaseWbegin = phaseWend = Phase::end;
    
		sur->addDmbrCom("Fields being maintained by specialist functions");
		sur->addDmbr(t_int,							"nverticesFrozen",	        "# of vertices on surface is frozen");
		sur->addDmbr(t_int,							"nvertices",		        "# of vertices on surface, change by calling MRISreallocVerticesAndFaces et al");
		sur->addDmbr(t_int,							"nfaces",		            "# of faces on surface, change by calling MRISreallocVerticesAndFaces et al");
		sur->addDmbr(t_bool,						"faceAttachmentDeferred",	"defer connecting faces to vertices for performance reasons");
		sur->addDmbr(t_int,							"nedges",		            "# of edges on surface");
		sur->addDmbr(t_int,							"nstrips");
		sur->addDmbr(t_PR_VERTEX_TOPOLOGY,			"vertices_topology");
		sur->addDmbr(t_PR_VERTEX,					"vertices");
		sur->addDmbr(t_ppVoid,						"dist_storage",			    "the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored");
		sur->addDmbr(t_ppVoid,						"dist_orig_storage",		"the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored");
		sur->addDmbr(t_int,							"tempsAssigned",	        "State of various temp fields that can be borrowed if not already in use");

	phaseRbegin = Phase::TopologyM; phaseWbegin = phaseWend = Phase::end;
		sur->addDmbr(t_PR_FACE,						"faces");
		sur->addDmbr(t_PR_MRI_EDGE,					"edges");
		sur->addDmbr(t_PR_FaceNormCacheEntry,		"faceNormCacheEntries");
		sur->addDmbr(t_PR_FaceNormDeferredEntry,	"faceNormDeferredEntries");

	phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;
		sur->addDmbr(t_PR_STRIP,					"strips");

	phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;
		sur->addDmbr(t_float,						"xctr");
		sur->addDmbr(t_float,						"yctr");
		sur->addDmbr(t_float,						"zctr");
		sur->addDmbr(t_float,						"xlo");
		sur->addDmbr(t_float,						"ylo");
		sur->addDmbr(t_float,						"zlo");
		sur->addDmbr(t_float,						"xhi");
		sur->addDmbr(t_float,						"yhi");
		sur->addDmbr(t_float,						"zhi");
		sur->addDmbr(t_float,						"x0", "center of spherical expansion");
		sur->addDmbr(t_float,						"y0");
		sur->addDmbr(t_float,						"z0");
		sur->addDmbr(t_PVERTEX,						"v_temporal_pole");			// WEIRD THAT THESE ARE POINTERS TO VERTICES
		sur->addDmbr(t_PVERTEX,						"v_frontal_pole");
		sur->addDmbr(t_PVERTEX,						"v_occipital_pole");
		sur->addDmbr(t_float,						"max_curv");
		sur->addDmbr(t_float,						"min_curv");
		sur->addDmbr(t_float,						"total_area");
		sur->addDmbr(t_double,						"avg_vertex_area");

	phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::XYZPositionConsequencesM;
		sur->addDmbr(t_double,						"avg_vertex_dist",			"set by MRIScomputeAvgInterVertexDist");
		sur->addDmbr(t_double,						"std_vertex_dist");
		sur->addDmbr(t_float,						"orig_area");
		sur->addDmbr(t_float,						"neg_area");
		sur->addDmbr(t_float,						"neg_orig_area",			"amount of original surface in folds");
		
	phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;
		sur->addDmbr(t_int,							"zeros");
		sur->addDmbr(t_int,							"hemisphere",			"which hemisphere");
		
	phaseRbegin = Phase::ExistenceM; phaseWbegin = Phase::end; phaseWend = Phase::end;
		sur->addDmbr(t_int,							"initialized");

		sur->addDmbrList("LIST_OF_MRIS_ELTS_3");

		sur->addDmbr(t_PLTA,						"lta");
		sur->addDmbr(t_PMATRIX,						"SRASToTalSRAS_");
		sur->addDmbr(t_PMATRIX,						"TalSRASToSRAS_");
		sur->addDmbr(t_int,							"free_transform");
		sur->addDmbr(t_double,						"radius",				"radius (if status==MRIS_SPHERE)");
		sur->addDmbr(t_float,						"a");
		sur->addDmbr(t_float,						"b");
		sur->addDmbr(t_float,						"c",					"ellipsoid parameters");

	phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::ExistenceM;
		sur->addDmbr(t_MRIS_fname_t,				"fname",				"file it was originally loaded from");

	phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;
		sur->addDmbr(t_float,						"Hmin",					"min mean curvature");
		sur->addDmbr(t_float,						"Hmax",					"max mean curvature");
		sur->addDmbr(t_float,						"Kmin",					"min Gaussian curvature");
		sur->addDmbr(t_float,						"Kmax",					"max Gaussian curvature");
		sur->addDmbr(t_double,						"Ktotal",				"total Gaussian curvature");

	phaseRbegin = phaseWbegin = Phase::ExistenceM; phaseWend = Phase::ExistenceM;
		sur->addDmbr(t_MRIS_Status,					"status",				"type of surface (e.g. sphere,"" plane)");
		sur->addDmbr(t_MRIS_Status,					"origxyz_status",		"type of surface (e.g. sphere, plane) that this origxyz were obtained from");
		sur->addDmbr(t_int,							"patch",				"if a patch of the surface");

	phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;
		sur->addDmbr(t_int,							"nlabels");
		sur->addDmbr(t_PMRIS_AREA_LABEL,			"labels",				"nlabels of these (may be null)");
		
	phaseRbegin = phaseWbegin = Phase::end; phaseWend = Phase::end;
		sur->addDmbr(t_char,						"nsize",				"size of neighborhoods or -1");
		sur->addDmbr(t_uchar,						"vtotalsMightBeTooBig", "MRISsampleDistances sets this");
		sur->addDmbr(t_short,						"nsizeMaxClock",		"changed whenever an edge is added or removed, which invalidates the vertex v#num values")->setNoHash();
		sur->addDmbr(t_char,						"max_nsize",			"max the neighborhood size has been set to (typically 3)");
		sur->addDmbr(t_char,						"dist_nsize",			"max mrisComputeVertexDistances has computed distances out to");
		sur->addDmbr(t_char,						"dist_orig_nsize",		"max mrisComputeOriginalVertexDistances has computed distances out to");
		sur->addDmbr(t_char,						"dist_alloced_flags",	"two flags, set when any dist(1) or dist_orig(2) allocated");
		sur->addDmbr(t_float,						"avg_nbrs",				"mean # of vertex neighbors");

	phaseRbegin = phaseWbegin = Phase::XYZPositionM; phaseWend = Phase::end;
		sur->addDmbr(t_pVoid,						"vp",					"for misc. use");
		sur->addDmbr(t_float,						"alpha",				"rotation around z-axis");
		sur->addDmbr(t_float,						"beta",					"rotation around y-axis");
		sur->addDmbr(t_float,						"gamma",				"rotation around x-axis");
		sur->addDmbr(t_float,						"da");
		sur->addDmbr(t_float,						"db");
		sur->addDmbr(t_float,						"dg",					"old deltas");
		sur->addDmbr(t_int,							"type",					"what type of surface was this initially");

	phaseRbegin = Phase::ExistenceM; phaseWbegin = Phase::end; phaseWend = Phase::end;
		sur->addDmbr(t_int,							"max_vertices",			"may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces");
		sur->addDmbr(t_int,							"max_faces",			"may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces");

		sur->addDmbr(t_MRIS_subject_name_t,			"subject_name",			"name of the subject");
		sur->addDmbr(t_float,						"canon_area");
		sur->addDmbr(t_int,							"noscale",				"don't scale by surface area if true");
		sur->addDmbr(t_PR_float,					"dx2",					"an extra set of gradient (not always alloced)");
		sur->addDmbr(t_PR_float,					"dy2");
		sur->addDmbr(t_PR_float,					"dz2");
		sur->addDmbr(t_PCOLOR_TABLE,				"ct");
		sur->addDmbr(t_int,							"useRealRAS",			"if 0 (default), vertex position is a conformed volume RAS with c_(r,\"a\",\"s\")=0.  "
																			"else is a real RAS (volume stored RAS)");
		sur->addDmbr(t_VOL_GEOM,					"vg",					"volume info from which this surface is created. valid iff vg.valid = 1");
		sur->addDmbr(t_MRIS_cmdlines_t,				"cmdlines")->setNoHash();
		sur->addDmbr(t_int,							"ncmds");
		sur->addDmbr(t_float,						"group_avg_surface_area",		"average of total surface area for group");
		sur->addDmbr(t_int,							"group_avg_vtxarea_loaded",		"average vertex area for group at each vertex");
		sur->addDmbr(t_int,							"triangle_links_removed",		"for quad surfaces");
		sur->addDmbr(t_pVoid,						"user_parms",					"for whatever the user wants to hang here");
		sur->addDmbr(t_PMATRIX,						"m_sras2vox",					"for converting surface ras to voxel");
		sur->addDmbr(t_PMRI,						"mri_sras2vox",					"volume that the above matrix is for");
		sur->addDmbr(t_pVoid,						"mht");
		sur->addDmbr(t_pVoid,						"temps");

		sur->addDmbrList("LIST_OF_MRIS_ELTS");
		sur->addDmbrListSublist("LIST_OF_MRIS_ELTS_1");
		sur->addDmbrListSublist("LIST_OF_MRIS_ELTS_3");
		
		Built built;
		built.push_back(fac);
		built.push_back(vtx);
		built.push_back(sur);
		return built;
	}

}
