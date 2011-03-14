# Makefile for proctrack/linuxproc plugin

AUTOMAKE_OPTIONS = foreign

PLUGIN_FLAGS = -module -avoid-version --export-dynamic 

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/src/common

pkglib_LTLIBRARIES = proctrack_cgroup.la

# Process group ID process tracking plugin.
proctrack_cgroup_la_SOURCES = proctrack_cgroup.c
proctrack_cgroup_la_LIBADD = $(top_builddir)/src/common/libcommon.la
proctrack_cgroup_la_LDFLAGS = $(SO_LDFLAGS) $(PLUGIN_FLAGS)
