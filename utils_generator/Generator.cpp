// Compile with     g++ -std=c++11 Generator.cpp
// Run with         rm -rf tmp ; mkdir tmp ; ./a.out ./tmp/
// Merge with       bcompare ./tmp ../include &
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
        string which;
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
        Dmbr* setWhich(string const & to) { which = to; return this; }
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
		}
		void os_push(ostream * new_tos) { os.push_back(new_tos); }
		void os_pop() { os.pop_back(); }
	private:
		vector<ostream *> os;
	};

	// Generate the various representations
	//
	struct IsImplementedSupplier {
		// Representations may only supply some of the properties
		virtual bool isImplemented(Clss const & c, Dmbr const & d) { return true; }
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
			indent() << "#define SEPARATE_VERTEX_TOPOLOGY" << endl;
			for (auto & cp : built) generateClass(*cp);
			for (auto & cp : built) generateMacros(*cp);
		}
	};


	struct Mrispv_implemented_properties : public IsImplementedSupplier {
		set<string> implemented_names;
		Mrispv_implemented_properties() {
			static const char* implemented_properties[] = {
				// in			for the metric properties calculations
				"status",
				"origxyz_status",
				"nvertices",    "vertices",
				"nfaces",       "faces",
				"nsize",
				"radius",
				// in out		for the metric properties calculations
				"dist_nsize",
				// out			for the metric properties calculations
				"xlo",
				"xhi",
				"ylo",
				"yhi",
				"zlo",
				"zhi",
				"xctr",
				"yctr",
				"zctr",
				"total_area",
				"avg_vertex_area",
				"avg_vertex_dist",
				"std_vertex_dist",
				"neg_orig_area",
				"neg_area",
				// in			for the metric properties calculations
				"v_ripflag",
                "v_num",
                "v_f",
				// in out		for the metric properties calculations
				"v_dist_capacity",
				"v_border",
				"v_x",
				"v_y",
				"v_z",
				"v_origarea",
				// out			for the metric properties calculations
				"v_area",
				"v_nx",
				"v_ny",
				"v_nz",
				"v_neg",
				"v_dist",
				// in			for the metric properties calculations 
				"f_ripflag",
                "f_v",
				"f_norm_orig_area",
				// in out		for the metric properties calculations
				// out			for the metric properties calculations
				"f_area",
				"f_normSet",
				"f_norm",
				"f_angle",
				nullptr };

			for (auto i = implemented_properties; *i; i++) {
				implemented_names.insert(*i);
			}

			// "v_VSize",
			// 	ELT(const, VERTEX_TOPOLOGY const *, vertices_topology) \
			// 	ELTX(const, FACE_TOPOLOGY   const *, faces_topology)
		}

		struct ImplementedProperty {
			bool   isImplemented, isVector;
			string name;
		};

		ImplementedProperty property(Clss const & c, Dmbr const & d) {
			string vectorPrefix;
			if (c.id == "Vertex")  vectorPrefix = "v_";		else
			if (c.id == "Face")    vectorPrefix = "f_";		else
			if (c.id == "Surface") vectorPrefix = "";		else
									assert(false);

			ImplementedProperty ip;
			ip.isImplemented = false;
			ip.isVector = vectorPrefix.size() > 0;
			ip.name = vectorPrefix + d.id;

			ip.isImplemented = (implemented_names.find(ip.name) != implemented_names.end());

			return ip;
		}

		virtual bool isImplemented(Clss const & c, Dmbr const & d) {
			return property(c, d).isImplemented;
		}
	} mrispv_implemented_properties;


	struct Generate_mrispv : public Generate_Base {

		void generateClass(Clss & c) {
			ColumnedOutput::T cols;

			for (auto & d : c.dmbr) {

				if (!d->type && d->commentNature != ComGeneral) continue;

				const char * comment = "//";
				if (d->type) {

					auto property = mrispv_implemented_properties.property(c, *d);
					if (!property.isImplemented) continue;

					cols << d->type->id << endC;
					if (property.isVector) cols << "*" << endC;
					cols << Just_left << property.name << endC << ";";
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
		}

		void generateMacros(Clss & c) {
		}

		Generate_mrispv(
			ostream & os,
			Built const & built) : Generate_Base(os, built)
		{

			indent() << "struct MRISPV {" << endl;
			depth++;
			for (auto & cp : built) generateClass(*cp);
			depth--;
			indent() << "};" << endl;

			for (auto & cp : built) generateMacros(*cp);
		}
	};

	// Generating the various Surface Vector Face accessors
	//
	struct Generate_abstract_spec_Base;
    struct Generate_abstract_inco_Base;

	struct Generate_SurfaceFromMRIS_inco;
	struct Generate_SurfaceFromMRIS_spec;
	struct Generate_SurfaceFromMRIS_impl;

	struct Generate_SurfaceFromMRISPV_inco;
	struct Generate_SurfaceFromMRISPV_spec;
	struct Generate_SurfaceFromMRISPV_impl;

	struct Generate_abstract_Base : public Generate_Base {

		std::string absolutePns(Phase::T p) {
			return "SurfaceFromMRIS_" + Phase::namespaceName(p);
		}

		virtual Generate_abstract_spec_Base*	 toGenerate_abstract_spec_Base    () { return nullptr; }
		virtual Generate_abstract_inco_Base*	 toGenerate_abstract_inco_Base    () { return nullptr; }

		virtual Generate_SurfaceFromMRIS_spec*   toGenerate_SurfaceFromMRIS_spec  () { return nullptr; }
		virtual Generate_SurfaceFromMRIS_impl*   toGenerate_SurfaceFromMRIS_impl  () { return nullptr; }
		virtual Generate_SurfaceFromMRIS_inco*   toGenerate_SurfaceFromMRIS_inco  () { return nullptr; }

		virtual Generate_SurfaceFromMRISPV_spec* toGenerate_SurfaceFromMRISPV_spec() { return nullptr; }
		virtual Generate_SurfaceFromMRISPV_impl* toGenerate_SurfaceFromMRISPV_impl() { return nullptr; }
		virtual Generate_SurfaceFromMRISPV_inco* toGenerate_SurfaceFromMRISPV_inco() { return nullptr; }

		virtual string classPrefix(Clss & c) = 0;
		virtual void   beginClass (Phase::T p, Clss & c) = 0;
		virtual void   endClass	  (Phase::T p, Clss & c) = 0;

		struct IsImplemented {
			IsImplementedSupplier* supplier;
			IsImplemented() : supplier(nullptr) {}
		} isImplementedSupplier;
		bool isImplemented(Clss const & c, Dmbr const & d) const {
			return isImplementedSupplier.supplier->isImplemented(c, d);
		}

		string baseClass(Phase::T p, Clss& c) {
#if 1
			return "Repr_Elt";
#else
			string result;
			auto const enhances = Phase::enhances[p];
			if (enhances == Phase::end) {
				result = "Repr_Elt";
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
						if (toGenerate_abstract_spec_Base()) cols << endC; 
						else if (needsSep) cols << " ";
					};

					auto maybeJust = [&](Justify j) {
						return toGenerate_abstract_spec_Base() ? j : Just_none;
					};
					
					vector<Dmbr*> deferredComGeneral;

					for (auto & dp : c.dmbr) {
						auto const & d = *dp;

						if (d.isGeneralComment()) {
							if (toGenerate_abstract_spec_Base()) {
								deferredComGeneral.push_back(dp);
							}
							continue;
						}

                        if (c.id == "Surface" && d.id == "vertices") bpt();

						if (!d.type) continue;

						if (!isImplemented(c, d)) continue;

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

						if (toGenerate_abstract_spec_Base()) 
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

                    if (!write && c.id == "Vertex") {
						if (toGenerate_abstract_spec_Base()) 
							cols << "inline ";
						else
							cols << ignoreWidth;
                        cols << "void";
                        maybeEndC();
                        cols << classPrefix(c) << maybeJust(Just_left) << "which_coords";
                        maybeEndC(false); 
                        cols << ignoreWidth
                             << "(int which, float *x, float *y, float *z) const";
                        maybeEndC();
                        generate_which_coords(p, c, cols, write);
                        cols << endR;
                    }                    
				} cols.append(tos(),depth);
			}
			endClass(p, c);
		}

        virtual void generateBeforeComment(ColumnedOutput::T& cols, bool const write) {
			cols << ";";
        }
        
		virtual void generateBeforeComment(ColumnedOutput::T& cols, Dmbr const & d, bool const write) {
		    generateBeforeComment(cols, write);
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, Clss const & c, Dmbr const & d, bool const write) {
		}

        virtual void generate_which_coords(Phase::T p, Clss & c, ColumnedOutput::T& cols, bool const write) {
			generateBeforeComment(cols, write);
        }

		virtual void generateNamespaceMacros(Phase::T p) {
		}

		void generateNamespace(string const & parent, string const & name, string const & representation) {
			indent() << "namespace " << name << " {" << endl;
			depth++;
			indent() << "typedef " << representation << " Representation;" << endl;

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
						generateNamespaceMacros(p);
						for (auto & cp : built) generateClass(p, *cp);
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
					for (auto & cp : built) {
						indent() << "friend struct " << name << "::" << Phase::namespaceName(p) << "::" << cp->id << ";" << endl;
					}
				}
				depth--;
				indent() << "};" << endl;
			}

			depth--;
			indent() << "} // " << name << endl;
		}

		Generate_abstract_Base(ostream & os, Built const & built) : Generate_Base(os, built)
		{
		}
	};

	struct Generate_abstract_inco_Base : public Generate_abstract_Base {
        virtual Generate_abstract_inco_Base*	 toGenerate_abstract_inco_Base    () { return this; }

		Generate_abstract_inco_Base(
			ostream & os,
            string const & parent,
			Built const & built) : Generate_abstract_Base(os, built)
		{
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

	struct Generate_abstract_spec_Base : public Generate_abstract_Base {

		virtual Generate_abstract_spec_Base* toGenerate_abstract_spec_Base() { return this; }

		Generate_abstract_spec_Base(
            ostream & os,
            string const & parent,
            Built const & built) : Generate_abstract_Base(os, built)
		{
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
			cols << "inline " << c.id << endC << "(" << endC << "Representation* representation" << (isSurface?"":", size_t idx") << endC << ");" << endR;

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

	struct Generate_abstract_impl_Base : public Generate_abstract_Base {

		Generate_abstract_impl_Base(
            ostream & os, 
            string const & parent,
            Built const & built) : Generate_abstract_Base(os, built) {
		}

		virtual string classPrefix(Clss & c) {
			return c.id + "::";
		}

		virtual void beginClass(Phase::T p, Clss & c) {
			ColumnedOutput::T cols;
            
            bool isSurface = (c.id == "Surface");

			cols << classPrefix(c) << c.id << endC << "(" << endC << "" << endC << ") {}" << endR;
			cols << classPrefix(c) << c.id << endC << "(" << endC << "Representation* representation" << (isSurface?"":", size_t idx") << endC 
                << ") : " << baseClass(p, c) << "(representation," << (isSurface  ? "0" : "idx" ) << ") {}" << endR;
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

	};


	// SurfaceFromMRIS
	//
	struct Generate_SurfaceFromMRIS_IsImplemented : public IsImplementedSupplier {
	} generate_SurfaceFromMRIS_IsImplemented;

	struct Generate_SurfaceFromMRIS_inco : public Generate_abstract_inco_Base {
		virtual Generate_SurfaceFromMRIS_inco* toGenerate_SurfaceFromMRIS_inco() { return this; }
		Generate_SurfaceFromMRIS_inco(
			ostream & os,
			string const & parent,
			Built const & built) : Generate_abstract_inco_Base(os, parent, built)
		{
			isImplementedSupplier.supplier = &generate_SurfaceFromMRIS_IsImplemented;
			generateNamespace(parent, "SurfaceFromMRIS", "MRIS");
		}
	};

	struct Generate_SurfaceFromMRIS_spec : public Generate_abstract_spec_Base {
		virtual Generate_SurfaceFromMRIS_spec* toGenerate_SurfaceFromMRIS_spec() { return this; }

		Generate_SurfaceFromMRIS_spec(
			ostream & os,
			string const & parent,
			Built const & built) : Generate_abstract_spec_Base(os, parent, built)
		{
			isImplementedSupplier.supplier = &generate_SurfaceFromMRIS_IsImplemented;
			generateNamespace(parent, "SurfaceFromMRIS", "MRIS");
		}

	};

	struct Generate_SurfaceFromMRIS_impl : public Generate_abstract_impl_Base {
		virtual Generate_SurfaceFromMRIS_impl* toGenerate_SurfaceFromMRIS_impl() { return this; }

		Generate_SurfaceFromMRIS_impl(
			ostream & os,
			string const & parent,
			Built const & built) : Generate_abstract_impl_Base(os, parent, built)
		{
			isImplementedSupplier.supplier = &generate_SurfaceFromMRIS_IsImplemented;
			generateNamespace(parent, "SurfaceFromMRIS", "MRIS");
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, Clss const & c, Dmbr const & d, bool const write) {
			auto t = d.type;
			auto tr = t->toPointerToRepeatedType();
			cols << ignoreWidth << Just_left << "    ";

			string raw = "repr->";
			if (d.mrisVector) raw += d.mrisVector, raw += "[idx].";
			raw += d.id;
			if (tr) raw += "[i]";

			if (c.id == "Face" && d.id == "v") {
				if (!write) {
					cols << "return Vertex(repr,repr->faces[idx].v[i])";
				}
				else {
					cols << "cheapAssert(repr == to.repr); repr->faces[idx].v[i] = to.idx";
				}
			}
			else if (c.id == "Vertex" && d.id == "f") {
				if (!write) {
					cols << "return Face(repr," << raw << ")";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << " = to.idx";
				}
			}
			else if (c.id == "Vertex" && d.id == "v") {
				if (!write) {
					cols << "return Vertex(repr," << raw << ")";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << " = to.idx";
				}
			}
			else if (c.id == "Surface" && (d.id == "v_temporal_pole" || d.id == "v_frontal_pole" || d.id == "v_occipital_pole")) {
				if (!write) {
					cols << "return Vertex(repr," << raw << " - repr->vertices)";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << " = repr->vertices + to.idx";
				}
			}
			else if (c.id == "Surface" && d.id == "vertices") {
                bpt();
				if (!write) {
					cols << "return Vertex(repr,i)";
				}
				else {
					assert(false);
				}
			}
			else if (c.id == "Surface" && d.id == "faces") {
				if (!write) {
					cols << "return Face(repr,i)";
				}
				else {
					assert(false);
				}
			}
			else {
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

        virtual void generate_which_coords(Phase::T p, Clss & c, ColumnedOutput::T& cols, bool const write) {
			cols << ignoreWidth << Just_left << "{" << endR;
            cols << ignoreWidth << Just_left << "    " << endC;
			cols << "MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);" << endR;
            cols << ignoreWidth << Just_left << "}" << endR;
        }
	};

	// SurfaceFromMRISPV
	//
	struct Generate_SurfaceFromMRISPV_inco : public Generate_abstract_inco_Base {
		virtual Generate_SurfaceFromMRISPV_inco* toGenerate_SurfaceFromMRISPV_inco() { return this; }
		Generate_SurfaceFromMRISPV_inco(
			ostream & os,
			string const & parent,
			Built const & built) : Generate_abstract_inco_Base(os, parent, built)
		{
			isImplementedSupplier.supplier = &mrispv_implemented_properties;
			generateNamespace(parent, "SurfaceFromMRISPV", "MRISPV");
		}
	};

	struct Generate_SurfaceFromMRISPV_spec : public Generate_abstract_spec_Base {
		virtual Generate_SurfaceFromMRISPV_spec* toGenerate_SurfaceFromMRISPV_spec() { return this; }

		Generate_SurfaceFromMRISPV_spec(
			ostream & os,
			string const & parent,
			Built const & built) : Generate_abstract_spec_Base(os, parent, built)
		{
			isImplementedSupplier.supplier = &mrispv_implemented_properties;
			generateNamespace(parent, "SurfaceFromMRISPV", "MRISPV");
		}

	};

	struct Generate_SurfaceFromMRISPV_impl : public Generate_abstract_impl_Base {
		virtual Generate_SurfaceFromMRISPV_impl* toGenerate_SurfaceFromMRISPV_impl() { return this; }

		Generate_SurfaceFromMRISPV_impl(
			ostream & os,
			string const & parent,
			Built const & built) : Generate_abstract_impl_Base(os, parent, built)
		{
			isImplementedSupplier.supplier = &mrispv_implemented_properties;
			generateNamespace(parent, "SurfaceFromMRISPV", "MRISPV");
		}

		virtual void generateMoreLines(ColumnedOutput::T& cols, Clss const & c, Dmbr const & d, bool const write) {
			auto implemented = mrispv_implemented_properties.property(c, d);

			auto t = d.type;
			auto tr = t->toPointerToRepeatedType();
			cols << ignoreWidth << Just_left << "    ";

			string raw = "repr->" + implemented.name;
			if (implemented.isVector) raw += "[idx]";

			if (tr) raw += "[i]";

			if (c.id == "Face" && d.id == "v") {
				if (!write) {
					cols << "return Vertex(repr," << raw << "[i]" << ")";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << "[i]" << " = to.idx";
				}
			}
			else if (c.id == "Vertex" && d.id == "f") {
				if (!write) {
					cols << "return Face(repr," << raw << ")";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << " = to.idx";
				}
			}
			else if (c.id == "Vertex" && d.id == "v") {
				if (!write) {
					cols << "return Vertex(repr," << raw << ")";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << " = to.idx";
				}
			}
			else if (c.id == "Surface" && (d.id == "v_temporal_pole" || d.id == "v_frontal_pole" || d.id == "v_occipital_pole")) {
				if (!write) {
					cols << "return Vertex(repr," << raw << " - repr->vertices)";
				}
				else {
					cols << "cheapAssert(repr == to.repr); " << raw << " = repr->vertices + to.idx";
				}
			}
			else if (c.id == "Surface" && d.id == "vertices") {
                bpt();
				if (!write) {
					cols << "return Vertex(repr,i)";
				}
				else {
					assert(false);
				}
			}
			else if (c.id == "Surface" && d.id == "faces") {
				if (!write) {
					cols << "return Face(repr,i)";
				}
				else {
					assert(false);
				}
			}
			else {
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
        
        virtual void generate_which_coords(Phase::T p, Clss & c, ColumnedOutput::T& cols, bool const write) {
			cols << ignoreWidth << Just_left << "{" << endR;
            cols << ignoreWidth << Just_left << "    " << endR;
			// cols << "repr->vertexCoord2XYZ(idx, which, x, y, z);" << endR;
            
cols << "#define CASE(WHICH, FIELD) \\" << endR;
cols << "  case WHICH: \\" << endR;
cols << "    *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \\" << endR;
cols << "    break;" << endR;
cols << "" << endR;
cols << "  switch (which) {" << endR;

            for (auto & dp : c.dmbr) {
                auto const & d = *dp;
                if (d.which == "") continue;
                if (!isImplemented(c, d)) continue;  
                bool canRead  = (d.firstReadingPhase <= p);
                if (!canRead) continue;
                assert(d.type->id == "float");
cols << ("    CASE("+d.which+","+d.id.substr(0,d.id.size()-1)+")") << endR;
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

	// Generate
	//
	void generate(Built const & built) {
		
		// Generate the various representations
		//
		if (true) {
			ofstream os(createFilesWithin + "mrisurf_FACE_VERTEX_MRIS_generated.h");
			Generate_mris(os, built);
		}
		if (true) {
			ofstream os(createFilesWithin + "mrisurf_MRIS_PropertiesInVectors.h");
			Generate_mrispv(os, built);
		}

		// Generate the Surface Vector Face accessors to the various representations
		//
		if (true) {
			string const root_fnm = "mrisurf_SurfaceFromMRIS_generated";
			ofstream os(createFilesWithin + root_fnm + ".h");
            os << "#pragma once" << endl;
            
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

		if (true) {
			string const root_fnm = "mrisurf_SurfaceFromMRISPV_generated";
			ofstream os(createFilesWithin + root_fnm + ".h");
            os << "#pragma once" << endl;

			auto fnm_inco = root_fnm + "_prefix";
			os << "#include \"./" << fnm_inco << ".h\"" << endl;
			{
				ofstream os_inco(createFilesWithin + fnm_inco + ".h");
				Generate_SurfaceFromMRISPV_inco(os_inco, fnm_inco, built);
			}

			{
				Generate_SurfaceFromMRISPV_spec(os, root_fnm, built);
			}

			auto fnm_impl = root_fnm + "_suffix";
			os << "#include \"./" << fnm_impl << ".h\"" << endl;
			{
				ofstream os_impl(createFilesWithin + fnm_impl + ".h");
				Generate_SurfaceFromMRISPV_impl(os_impl, fnm_impl, built);
			}
		}

		auto cmd = "dir " + createFilesWithin;
		cout << cmd << endl;
		system(cmd.c_str());
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
		vtx->addDmbr(t_float,		"x"					, "current coordinates	")->setWhich("CURRENT_VERTICES");
		vtx->addDmbr(t_float,		"y"					, "use MRISsetXYZ() to set");
		vtx->addDmbr(t_float,		"z"					);

		phaseRbegin = Phase::XYZPositionM; phaseWbegin = Phase::end; phaseWend = Phase::end;

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"origx"				, "original coordinates, see also MRIS::origxyz_status")->setWhich("ORIGINAL_VERTICES");
		vtx->addDmbr(t_float,		"origy"				, "use MRISsetOriginalXYZ(, ");
		vtx->addDmbr(t_float,		"origz"				, "or MRISsetOriginalXYZfromXYZ to set");

		phaseRbegin = phaseWbegin = Phase::XYZPositionConsequencesM; phaseWend = Phase::end;

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"nx")->setWhich("VERTEX_NORMALS");
		vtx->addDmbr(t_float,		"ny");
		vtx->addDmbr(t_float,		"nz", "curr normal");

		phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		vtx->addDmbr(t_float,		"pnx")->setWhich("PIAL_NORMALS");
		vtx->addDmbr(t_float,		"pny");
		vtx->addDmbr(t_float,		"pnz", "pial normal");

		vtx->addDmbrCom("");
		vtx->addDmbr(t_float,		"wnx")->setWhich("WHITE_NORMALS"); 
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

		vtx->addDmbr(t_float,		"cx")->setWhich("CANONICAL_VERTICES");
		vtx->addDmbr(t_float,		"cy"); 
		vtx->addDmbr(t_float,		"cz", "coordinates in canonical coordinate system");

		phaseRbegin = phaseWbegin = Phase::DistortM; phaseWend = Phase::end;

		vtx->addDmbr(t_float,		"tx")->setWhich("TMP_VERTICES");
		vtx->addDmbr(t_float,		"ty"); 
		vtx->addDmbr(t_float,		"tz", "tmp coordinate storage");
		vtx->addDmbr(t_float,		"t2x")->setWhich("TMP2_VERTICES"); 
		vtx->addDmbr(t_float,		"t2y"); 
		vtx->addDmbr(t_float,		"t2z", "another tmp coordinate storage");
		vtx->addDmbr(t_float,		"targx"); 
		vtx->addDmbr(t_float,		"targy"); 
		vtx->addDmbr(t_float,		"targz", "target coordinates");
		vtx->addDmbr(t_float,		"pialx")->setWhich("PIAL_VERTICES"); 
		vtx->addDmbr(t_float,		"pialy"); 
		vtx->addDmbr(t_float,		"pialz", "pial surface coordinates");
		vtx->addDmbr(t_float,		"whitex")->setWhich("WHITE_VERTICES"); 
		vtx->addDmbr(t_float,		"whitey"); 
		vtx->addDmbr(t_float,		"whitez", "white surface coordinates");
		vtx->addDmbr(t_float,		"l4x"); 
		vtx->addDmbr(t_float,		"l4y"); 
		vtx->addDmbr(t_float,		"l4z", "layerIV surface coordinates");
		vtx->addDmbr(t_float,		"infx")->setWhich("INFLATED_VERTICES"); 
		vtx->addDmbr(t_float,		"infy"); 
		vtx->addDmbr(t_float,		"infz", "inflated coordinates");
		vtx->addDmbr(t_float,		"fx")->setWhich("FLATTENED_VERTICES"); 
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
		sur->addDmbr(t_ppVoid,						"dist_storage",			    "the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored")->setNoHash();
		sur->addDmbr(t_ppVoid,						"dist_orig_storage",		"the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored")->setNoHash();
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
		sur->addDmbr(t_pVoid,						"vp",					"for misc. use")->setNoHash();
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
		sur->addDmbr(t_pVoid,						"user_parms",					"for whatever the user wants to hang here")->setNoHash();
		sur->addDmbr(t_PMATRIX,						"m_sras2vox",					"for converting surface ras to voxel");
		sur->addDmbr(t_PMRI,						"mri_sras2vox",					"volume that the above matrix is for");
		sur->addDmbr(t_pVoid,						"mht")->setNoHash();
		sur->addDmbr(t_pVoid,						"temps")->setNoHash();

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
