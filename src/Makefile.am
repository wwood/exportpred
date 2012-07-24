LDADD = ../lib/libghmm.la

bin_PROGRAMS=test_ghmm exportpred simulate_signalseqs

test_ghmm_CXXFLAGS = @CXXFLAGS@ @PCRE_CFLAGS@
test_ghmm_LIBS = @LIBS@ @PCRE_LIBS@
test_ghmm_SOURCES = test_ghmm.cc

exportpred_CXXFLAGS = @CXXFLAGS@ @PCRE_CFLAGS@
exportpred_LIBS = @LIBS@ @PCRE_LIBS@
exportpred_SOURCES = predict_pexel.cc predict_pexel.hh ss_model.cc signalp_model.cc

simulate_signalseqs_CXXFLAGS = @CXXFLAGS@ @PCRE_CFLAGS@
simulate_signalseqs_LIBS = @LIBS@ @PCRE_LIBS@
simulate_signalseqs_SOURCES = simulate_signalseqs.cc ss_model.cc signalp_model.cc
