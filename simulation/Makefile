# ============================================
#  Makefile para eo-silo con Box2D 2.4.2
#  Estructura: src/, obj/, bin/
# ============================================

# --- Nombres y Directorios ---
TARGET_NAME = eo-silo
SRC_DIR     = src
OBJ_DIR     = obj
BIN_DIR     = bin

# Ruta relativa a Box2D (ajustar si es necesario)
BOX2D       = ../../box2d

# --- Compilador y Flags ---
CXX      = g++
INCLUDES = -I$(BOX2D)/include -I$(SRC_DIR)
LIBPATH  = -L$(BOX2D)/build/bin

# Flags de compilación
CXXFLAGS = -O3 -std=gnu++17 -Wall $(INCLUDES)

# Flags de linkeo
LDFLAGS  = -O3 -std=gnu++17 -Wall $(LIBPATH) -no-pie
LIBS     = -lbox2d -lgsl -lgslcblas -lm

# --- Archivos ---
# Busca automáticamente todos los .cpp dentro de la carpeta src
SRCS = $(filter-out $(SRC_DIR)/save_tensors.cpp, $(wildcard $(SRC_DIR)/*.cpp))
# save_tensors.cpp: fragmento suelto (sin includes) que no llama nadie; no compila solo

# Convierte la lista de sources (src/%.cpp) en objetos (obj/%.o)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

# Ruta final del ejecutable
TARGET = $(BIN_DIR)/$(TARGET_NAME)

# ============================================
#  Reglas
# ============================================

.PHONY: all clean directories

all: directories $(TARGET)

# Regla para crear directorios si no existen
directories:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

# Linkeo final
$(TARGET): $(OBJS)
	@echo "Linkeando $(TARGET)..."
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
	@echo "Construcción exitosa. Ejecutable en: $(TARGET)"

# Compilación de objetos (de src/ a obj/)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compilando $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Limpieza
clean:
	@echo "Limpiando archivos compilados..."
	rm -rf $(OBJ_DIR) $(BIN_DIR)