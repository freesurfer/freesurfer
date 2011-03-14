# Makerules.
# This file is part of AutoTroll.
# Copyright (C) 2006, 2007, 2009, 2010  Benoit Sigoure.
#
# AutoTroll is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
# USA.
#
# In addition, as a special exception, the copyright holders of AutoTroll
# give you unlimited permission to copy, distribute and modify the configure
# scripts that are the output of Autoconf when processing the macros of
# AutoTroll.  You need not follow the terms of the GNU General Public License
# when using or distributing such scripts, even though portions of the text of
# AutoTroll appear in them. The GNU General Public License (GPL) does govern
# all other use of the material that constitutes AutoTroll.
#
# This special exception to the GPL applies to versions of AutoTroll
# released by the copyright holders of AutoTroll.  Note that people who make
# modified versions of AutoTroll are not obligated to grant this special
# exception for their modified versions; it is their choice whether to do so.
# The GNU General Public License gives permission to release a modified version
# without this exception; this exception also makes it possible to release a
# modified version which carries forward this exception.

 # ------------- #
 # DOCUMENTATION #
 # ------------- #

# See autotroll.m4 :)


SUFFIXES = .moc.cpp .moc.cc .moc.cxx .moc.C .h .hh .hpp \
           .ui .ui.h .ui.hh .ui.hpp \
           .qrc .qrc.cpp .qrc.cc .qrc.cxx .qrc.C

# --- #
# MOC #
# --- #

.hpp.moc.cpp:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.hh.moc.cpp:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.h.moc.cpp:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@

.hpp.moc.cc:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.hh.moc.cc:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.h.moc.cc:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@

.hpp.moc.cxx:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.hh.moc.cxx:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.h.moc.cxx:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@

.hpp.moc.C:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.hh.moc.C:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
.h.moc.C:
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@

%.moc: %.cpp
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@
	
moc_%.cpp: %.h
	$(MOC) $(QT_CPPFLAGS) $(EXTRA_CPPFLAGS) $< -o $@


# --- #
# UIC #
# --- #

.ui.ui.hpp:
	$(UIC) $< -o $@

.ui.ui.hh:
	$(UIC) $< -o $@

.ui.ui.h:
	$(UIC) $< -o $@

ui_%.h: %.ui
	$(UIC) $< -p -o $@


# --- #
# RCC #
# --- #

.qrc.qrc.cpp:
	$(RCC) -name `echo "$<" | sed 's|^.*/\(.*\)\.qrc$$|\1|'` $< -o $@

.qrc.qrc.cc:
	$(RCC) -name `echo "$<" | sed 's|^.*/\(.*\)\.qrc$$|\1|'` $< -o $@

.qrc.qrc.cxx:
	$(RCC) -name `echo "$<" | sed 's|^.*/\(.*\)\.qrc$$|\1|'` $< -o $@

.qrc.qrc.C:
	$(RCC) -name `echo "$<" | sed 's|^.*/\(.*\)\.qrc$$|\1|'` $< -o $@

qrc_%.cpp: %.qrc
	$(RCC) $<  -o $@

DISTCLEANFILES = $(BUILT_SOURCES)
