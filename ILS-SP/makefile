CPLEX_VERSION = 201

#detecta se o sistema é de 32 ou 64 bits
BITS_OPTION = -m64

### diretorios do boost
BOOSTDIR = /home/carlos/boost_1_74_0

####diretorios com as libs do cplex
CPLEXDIR  = /opt/ibm/ILOG/CPLEX_Studio$(CPLEX_VERSION)/cplex
CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio$(CPLEX_VERSION)/concert
   
CPLEXLIBDIR   = $(CPLEXDIR)/lib/x86-64_linux/static_pic
CONCERTLIBDIR = $(CONCERTDIR)/lib/x86-64_linux/static_pic

#### define o compilador
CPPC = g++
#############################

#### opcoes de compilacao e includes
#CCOPT = $(BITS_OPTION) -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++2a -DROBUST_TWx -DROBUST_DEMANDx -DTW -DTESTINGx -DUSEMMDx -DUSEMSDx -DUSEFSRx  -g
CCOPT = $(BITS_OPTION) -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++2a -DROBUST_TWx -DROBUST_DEMANDx -DTWx -DTESTINGx -DUSEMMDx -DUSEMSDx -DUSEFSRx  -DDIMACSx -DSP_STATS -g
# CCOPT = $(BITS_OPTION) -O3 -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++17 -g -DROBUST_TW
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)
#############################

#### flags do linker
#CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -lconcert -lilocplex -lcplex -lm -lpthread -ldl
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -lpthread -lconcert -lilocplex -lcplex -lm -lpthread -ldl
#############################

#### diretorios com os source files e com os objs files
SRCDIR = src
OBJDIR = obj
#############################

#### lista de todos os srcs e todos os objs
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(SRCS))
#############################

#### regra principal, gera o executavel
ilssp.out: $(OBJS) 
	@echo  "\033[31m \nLinking all objects files: \033[0m"
	$(CPPC) $(BITS_OPTION) $(OBJS) -o $@ $(CCLNFLAGS)
############################

#inclui os arquivos de dependencias
-include $(OBJS:.o=.d)

#regra para cada arquivo objeto: compila e gera o arquivo de dependencias do arquivo objeto
#cada arquivo objeto depende do .c e dos headers (informacao dos header esta no arquivo de dependencias gerado pelo compiler)
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo  "\033[31m \nCompiling $<: \033[0m"
	$(CPPC) $(CCFLAGS) -c $< -o $@
	@echo  "\033[32m \ncreating $< dependency file: \033[0m"
	$(CPPC) -std=c++0x  -MM $< > $(basename $@).d
	@mv -f $(basename $@).d $(basename $@).d.tmp #proximas tres linhas colocam o diretorio no arquivo de dependencias (g++ nao coloca, surprisingly!)
	@sed -e 's|.*:|$(basename $@).o:|' < $(basename $@).d.tmp > $(basename $@).d
	@rm -f $(basename $@).d.tmp

#delete objetos e arquivos de dependencia
clean:
	@echo "\033[31mcleaning obj directory \033[0m"
	@rm ilssp.out -f $(OBJDIR)/*.o $(OBJDIR)/*.d


rebuild: clean ilssp.out
