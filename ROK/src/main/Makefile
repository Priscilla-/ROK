# Fundamental Variables
NAME = rok

# Extra Variables Required Independent of the System
OBJD = $(BASE)/obj/main
SRCD = $(BASE)/src/main
SRCD_KERR = $(BASE)/src/kerr
SRCD_AAF = $(BASE)/src/kerr/AAFormulation
SRCD_RESONACES = $(BASE)/src/resonances
SRCD_PNSF = $(BASE)/src/PNForce
SRCD_Ell = $(BASE)/src/elliptic_functions

INCLUDE = $(BASE)/include
EXEC = $(EXED)/$(NAME)

vpath

vpath %.o $(OBJD)
vpath %.c $(SRCD) $(SRCD_KERR) 
vpath %.c $(SRCD_Ell) $(SRCD_RESONACES) $(SRCD_PNSF) $(SRCD_AAF) $(SRCD_PNSF)
vpath %.h $(INCLUDE)

OBJS_KERR =bulirsch_stoer.o\
  check_if_orbit_plunges.o\
	checks_errors.o\
  compute_Keplerian_period.o\
  compute_new_constants_of_motion.o\
  compute_ELCQ_dot.o\
  christoffel_symbols.o\
	ELCQ_to_peiz.o\
  evolution_RR.o\
  flat_space_XVA.o\
	geodesic_evolution.o\
  geodesic_RHS.o\
	initialization.o\
	pez_to_ELCQ.o\
  turning_points.o\
  waveform.o\
	wrap_up.o
    
OBJS_AAF = analytic_r_costh.o\
    analytic_derivatives.o\
    bulirsch_stoer_AAV.o\
    evolution_AA.o\
    evolution_step.o\
    geodesic_evolution_AV.o\
    AA_RHS.o\
    turning_points_derivatives.o\
    
OBJS_Resonances = compute_frequencies.o\
    find_resonances.o\
    search_resonant_frequencies.o\
    
OBJS_PNSF = dchi_dpsi_dphi_dt.o\
    PN_ELzQdot.o\
    
OBJS_main = algebra.o\
    checks_main_errors.o\
    conversion_factors.o\
    create_files.o\
    main.o\
    modify_computational_parameters.o\
    reading_parameters.o\
    save_run_parameters.o

OBJS_ell = EllipticIntegrals.o\
    JacobianEllipticFunctions.o
                
							
$(NAME): $(OBJS_KERR) $(OBJS_main) $(OBJS_ell) $(OBJS_Resonances) $(OBJS_AAF) $(OBJS_PNSF)
	@echo 
	@echo ------------------------------------------------------------------
	@echo  Building $(NAME_GAME) Executable Version $(VERSION)
	@echo  Name of the Executable: $(NAME)
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(OBJS_main) $(OBJS_KERR) $(OBJS_ell) $(OBJS_Resonances) $(OBJS_AAF) $(OBJS_PNSF) -o $(EXEC) $(LIBRARIES) $(LINK_LIBRARIES)
	
$(OBJS_KERR): %.o: %.c
	@echo ------------------------------------------------------------------
	@echo  [$(NAME_GAME) $(VERSION)] Building Object: $*.o
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(CFLAGS) $(SRCD_KERR)/$(*).c
    
$(OBJS_AAF): %.o: %.c
	@echo ------------------------------------------------------------------
	@echo  [$(NAME_GAME) $(VERSION)] Building Object: $*.o
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(CFLAGS) $(SRCD_AAF)/$(*).c
    
$(OBJS_Resonances): %.o: %.c
	@echo ------------------------------------------------------------------
	@echo  [$(NAME_GAME) $(VERSION)] Building Object: $*.o
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(CFLAGS) $(SRCD_RESONACES)/$(*).c
    
$(OBJS_PNSF): %.o: %.c
	@echo ------------------------------------------------------------------
	@echo  [$(NAME_GAME) $(VERSION)] Building Object: $*.o
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(CFLAGS) $(SRCD_PNSF)/$(*).c
    	    
$(OBJS_main): %.o: %.c
	@echo ------------------------------------------------------------------
	@echo  [$(NAME_GAME) $(VERSION)] Building Object: $*.o
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(CFLAGS) $(SRCD)/$(*).c
    
$(OBJS_ell): %.o: %.c
	@echo ------------------------------------------------------------------
	@echo  [$(NAME_GAME) $(VERSION)] Building Object: $*.o
	@echo ------------------------------------------------------------------
	cd $(OBJD); $(CC) $(CFLAGS) $(SRCD_Ell)/$(*).c
    
clean:
	cd $(OBJD); $(RM) *.o
	cd $(SRCD); $(RM) *~
	cd $(SRCD_KERR); $(RM) *~
	cd $(SRCD_AAF); $(RM) *~
	cd $(SRCD_PNSF); $(RM) *~
	cd $(SRCD_RESONACES); $(RM) *~
	cd $(SRCD_Ell); $(RM) *~
	$(RM) $(EXEC)
		