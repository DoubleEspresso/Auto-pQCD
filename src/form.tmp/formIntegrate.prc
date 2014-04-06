#procedure formIntegrate(expr)

  *// simple wrapper to python routine which collects all F,J integrals
  *// and makes finite.out, as well as divergent.out for this diagram.  

  #write <integration/tmp.out> "%E" expr
*  #system cd integration && python collect_integrands.py