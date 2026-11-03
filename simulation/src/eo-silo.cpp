/** \mainpage ts-2
 * \brief Simulación de descarga de un silo con un extremo del orificio elevado.
 * 
 * Programa que implementa la biblioteca Box2D para simular la
 * descarga de un silo con orificio oblicuo conteniendo material 
 * granular, a diferencia de la versión previa
 * tilted-silo.
 *
 * \author Manuel Carlevaro <manuel@iflysib.unlp.edu.ar>
 *
 * \version 1.0 
 *
 * \date 2024.08.05
 *
 * */


#include <iostream>
using namespace std;
#include <box2d/box2d.h>
#include "globalSetup.hpp"
#include <string>
#include <cmath>
#include <vector>
using std::sin; using std::cos; using std::fabs; using std::sqrt;
#include <iomanip>
#include <sstream>
#include <fstream>   // <-- NUEVO
#include <cstdlib>
using std::exit;
#include "siloAux.hpp"
#include "jam_output.hpp"
#include "random.hpp"
#include <queue>
#include <unordered_map>
#include <algorithm>



#define PI 3.1415926

// Estructura para guardar qué zonas están libres
struct ReinjectionZones {
    std::vector<int> validBins; // Índices de franjas Y con poca densidad
    float binHeight;
    float minY;
    bool computed; // Para saber si ya se calculó en este frame
};

// Función que llena la estructura (Solo se llama UNA VEZ por frame si es necesario)
static void computeDensityMap(b2World& world, 
                              ReinjectionZones& zones,
                              float minX, float maxX, 
                              float minY, float maxY) 
{
    const int numBins = 20; 
    int counts[numBins] = {0};
    float zoneHeight = maxY - minY;
    zones.binHeight = zoneHeight / (float)numBins;
    zones.minY = minY;
    zones.validBins.clear();
    
    int totalInZone = 0;

    // Escaneo O(N) una sola vez
    for (b2Body* b = world.GetBodyList(); b; b = b->GetNext()) {
        uintptr_t ud = b->GetUserData().pointer;
        if (ud == 0) continue;
        BodyData* inf = reinterpret_cast<BodyData*>(ud);
        if (!inf || !inf->isGrain || !inf->isIn) continue;

        b2Vec2 p = b->GetPosition();
        if (p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY) {
            int bin = (int)((p.y - minY) / zoneHeight * numBins);
            if (bin >= 0 && bin < numBins) {
                counts[bin]++;
                totalInZone++;
            }
        }
    }

    // Determinar bins válidos (< 5%)
    float threshold = 0.05f;
    if (totalInZone == 0) {
        for(int i=0; i<numBins; i++) zones.validBins.push_back(i);
    } else {
        // Buscar zonas vacías
        for (int i = 0; i < numBins; ++i) {
            if ((float)counts[i] / (float)totalInZone < threshold) {
                zones.validBins.push_back(i);
            }
        }
        // Fallback: si todo está lleno, buscar el mínimo
        if (zones.validBins.empty()) {
            int minVal = 1000000;
            for(int i=0; i<numBins; i++) if(counts[i] < minVal) minVal = counts[i];
            for(int i=0; i<numBins; i++) if(counts[i] == minVal) zones.validBins.push_back(i);
        }
    }
    zones.computed = true;
}


// Rompe el arco buscando, para cada "rayo" vertical, el grano más bajo
// (menor y) dentro de una pequeña ventana en x encima del orificio,
// y reinyecta solo esos granos. Es equivalente a un raycast vertical
// pero hecho geométricamente, sin depender del RayCast de Box2D.
static void unjamArch(b2World& world,
                      GlobalSetup* globalSetup,
                      RandomGenerator& rng,
                      float minX, float maxX, // Nuevos parámetros X
                      float minY, float maxY) // Nuevos parámetros Y
{
    // --- Geometría del orificio ---
    double theta_rad = globalSetup->caja.theta * PI / 180.0;
    float rHole = (float)globalSetup->caja.r;

    // Proyección horizontal del orificio (para theta=0 es ±rHole)
    float xMin = (float)(-rHole * std::cos(theta_rad));
    float xMax = (float)( rHole * std::cos(theta_rad));

    float yHole = (float)globalSetup->caja.y_max_hole;
    double R0   = globalSetup->granos[0]->radio;   // radio típico de grano

    // Piso a partir del cual disparamos el pseudo-rayo (sin límite superior)
    float yBottom = yHole - (float)(2.0 * R0);     // un poco por debajo del orificio
    if (yBottom < 0.0f) yBottom = 0.0f;            // nunca bajo cero

    // "Ancho" en x de cada pseudo-rayo (ventanita alrededor de x_col)
    float halfDx = (float)(0.6 * R0); // ~60% del radio

    // Número de columnas (pseudo-rayos)
    int numCols = 80;

    // Para evitar duplicados
    std::vector<b2Body*> toReinject;
    toReinject.reserve(numCols);

    auto addUnique = [&](b2Body* b) {
        for (b2Body* bb : toReinject) {
            if (bb == b) return;
        }
        toReinject.push_back(b);
    };

    // --- Recorremos columnas sobre el orificio ---
    for (int i = 0; i < numCols; ++i) {
        float alpha = (i + 0.5f) / (float)numCols;
        float xCol  = xMin + alpha * (xMax - xMin);

        b2Body* bestBody = nullptr;
        float bestY = 1e9f;

        // Buscar el grano más bajo dentro de la ventanita en x
        for (b2Body* b = world.GetBodyList(); b; b = b->GetNext()) {
            uintptr_t ud = b->GetUserData().pointer;
            if (ud == 0) continue;
            BodyData* inf = reinterpret_cast<BodyData*>(ud);
            if (!inf || !inf->isGrain || !inf->isIn) continue;

            b2Vec2 p = b->GetPosition();

            // Debe estar en la franja vertical alrededor del orificio
            if (p.y < yBottom) continue;  // sin límite superior: buscamos hasta encontrar el primero

            // Debe estar dentro de la columna (ventana en x)
            if (p.x < xCol - halfDx || p.x > xCol + halfDx) continue;

            // Nos quedamos con el más cercano al orificio (menor y)
            if (p.y < bestY) {
                bestY   = p.y;
                bestBody = b;
            }
        }

        if (bestBody) {
            addUnique(bestBody);
        }
    }

    // --- Reinyectar todos los granos seleccionados ---
    int numReinjected = 0;

    for (b2Body* body : toReinject) {
        if (!body) continue;

        uintptr_t ud = body->GetUserData().pointer;
        if (ud == 0) continue;
        BodyData* inf = reinterpret_cast<BodyData*>(ud);
        if (!inf || !inf->isGrain) continue;

        // OBTENER RADIO para no chocar con paredes
        b2Fixture* f = body->GetFixtureList();
        b2Shape* bs  = (b2Shape*) f->GetShape();
        float r      = bs->m_radius;

        // --- CAMBIO AQUÍ: Usar las coordenadas pasadas por parámetro ---
        // Usamos minX/maxX ajustados por el radio del grano actual
        double x = rng.getAB(minX + r, maxX - r); 
        double y = rng.getAB(minY, maxY);         

        body->SetTransform(b2Vec2((float)x, (float)y), 0.0f);
        body->SetLinearVelocity(b2Vec2(0.0f, 0.0f));
        body->SetAngularVelocity(0.0f);
        body->SetAwake(true);
        inf->isIn = true;

        ++numReinjected;
    }

    if (numReinjected > 0) {
        std::cout << "# Unjam pseudo-rayos: reinyectados "
                  << numReinjected
                  << " granos en Y=[" << minY << ", " << maxY << "]" // Debug info actualizada
                  << std::endl;
    } else {
        std::cout << "# Unjam pseudo-rayos: no se encontraron granos en el arco."
                  << std::endl;
    }
}

void saveAtascoFrame(b2World* world, GlobalSetup* globalSetup, int jamIndex) {
    std::stringstream ss;
    ss << "frames_" << globalSetup->dirID << "/atasco_" 
       << std::setfill('0') << std::setw(5) << jamIndex << ".xy";
    
    std::string fileName = ss.str();
    std::ofstream f(fileName.c_str());

    if (!f.is_open()) {
        std::cout << "ERROR: No se pudo crear el archivo de atasco." << std::endl;
        return;
    }

    f << std::scientific << std::setprecision(6);

    for (b2Body* b = world->GetBodyList(); b; b = b->GetNext()) {
        uintptr_t ud = b->GetUserData().pointer;
        // Nota: El silo estático a veces tiene userData, a veces no, 
        // pero en tu main le pusiste: siloD->gID = -100.
        // Si el puntero es nulo, lo saltamos, pero el silo debería tener datos.
        if (ud == 0) continue;
        
        BodyData* inf = reinterpret_cast<BodyData*>(ud);

        // --- CAMBIO IMPORTANTE ---
        // Antes filtrábamos: if (!inf->isGrain) continue;
        // Ahora permitimos que pase si es grano (isIn) O si es el silo/tapa (gID < 0)
        
        bool esGranoActivo = (inf->isGrain && inf->isIn);
        bool esEstructura = (inf->gID < 0); // Silo (-100) o Tapa (-110)

        if (!esGranoActivo && !esEstructura) continue;

        b2Fixture* fixture = b->GetFixtureList();
        // A veces un cuerpo tiene varios fixtures, aquí asumimos 1 principal para simplificar,
        // o iteramos si el silo está compuesto. Tu silo es 1 chain, así que ok.
        while (fixture) {
            b2Shape::Type shapeType = fixture->GetType();

            // 1. POLÍGONOS (Granos cuadrados)
            if (shapeType == b2Shape::e_polygon) {
                b2PolygonShape* poly = (b2PolygonShape*)fixture->GetShape();
                int count = poly->m_count;
                f << inf->gID << " " << count << " ";
                for (int i = 0; i < count; ++i) {
                    b2Vec2 wv = b->GetWorldPoint(poly->m_vertices[i]);
                    f << wv.x << " " << wv.y << " ";
                }
                f << inf->tipo << "\n";
            }
            // 2. CADENAS (El Silo es un ChainShape)
            else if (shapeType == b2Shape::e_chain) {
                b2ChainShape* chain = (b2ChainShape*)fixture->GetShape();
                int count = chain->m_count;
                // El renderizador dibuja líneas conectadas si nVerts > 2
                f << inf->gID << " " << count << " ";
                
                // Box2D almacena vértices locales en el chain, pero como el body
                // del silo suele estar en (0,0) sin rotación, local = world.
                // Aun así, usamos GetWorldPoint para ser consistentes.
                for (int i = 0; i < count; ++i) {
                    // Accedemos a los vértices usando la API de chain
                    b2Vec2 v = chain->m_vertices[i];
                    b2Vec2 wv = b->GetWorldPoint(v);
                    f << wv.x << " " << wv.y << " ";
                }
                // Tipo 0 para estructura (o el que use tu render para negro)
                f << "0\n"; 
            }
            // 3. BORDES / EDGE (La Tapa del piso es un EdgeShape)
            else if (shapeType == b2Shape::e_edge) {
                b2EdgeShape* edge = (b2EdgeShape*)fixture->GetShape();
                // Un edge siempre tiene 2 vértices
                f << inf->gID << " 2 ";
                
                b2Vec2 v1 = b->GetWorldPoint(edge->m_vertex1);
                b2Vec2 v2 = b->GetWorldPoint(edge->m_vertex2);
                
                f << v1.x << " " << v1.y << " " 
                  << v2.x << " " << v2.y << " " 
                  << "0\n"; // Tipo 0
            }
            // 4. CÍRCULOS (Granos redondos)
            else if (shapeType == b2Shape::e_circle) {
                b2CircleShape* circle = (b2CircleShape*)fixture->GetShape();
                b2Vec2 pos = b->GetWorldPoint(circle->m_p);
                f << inf->gID << " 1 " 
                  << pos.x << " " << pos.y << " " 
                  << circle->m_radius << " " 
                  << inf->tipo << "\n";
            }

            fixture = fixture->GetNext();
        }
    }
    
    f.close();
    std::cout << "# Snapshot de atasco (con silo) guardado: " << fileName << std::endl;
}


int main(int argc, char *argv[])
{
    cout << "# eo-silo" << endl;
    cout << "# v1.0 [2024.08.05]" << endl;
    string inputParFile(argv[1]);
    GlobalSetup *globalSetup = new GlobalSetup(inputParFile);
    INFO("gs OK")
    if (globalSetup->caja.theta < 0.0) {
        cout << "ERROR: en eo-silo solo se puede usar theta >= 0." << endl;
        exit(1);
    }
    RandomGenerator rng(globalSetup->randomSeed);
    string folder_cmd = "mkdir -p frames_"+globalSetup->dirID;
    int sysret = system(folder_cmd.c_str());
    cout << "# Creación de directorio de frames: " << (sysret ? "NOK." : "OK.") 
        << endl;
    cout << "# Creación del sistema ..." << endl;



    // Definición del mundo
    double theta_rad = globalSetup->caja.theta * PI / 180.0;
    b2Vec2 gravedad(0.0f, globalSetup->g);
    b2World world(gravedad); // Llenado con el silo vertical
    cout << "#\t- Objeto mundo creado." << endl;
    
    // Definición del contenedor
    b2BodyDef bd;
    bd.position.Set(0.0f, 0.0f);
    bd.type = b2_staticBody;
    BodyData *siloD = new BodyData;
    siloD->isGrain = false;
    siloD->nLados = 8;
    siloD->gID = -100;
    bd.userData.pointer = uintptr_t(siloD);
    b2Body *silo = world.CreateBody(&bd);
    b2Vec2 wall_poly[6];
    wall_poly[0].Set(-globalSetup->caja.r * cos(theta_rad), 0.0f);
    wall_poly[1].Set(-globalSetup->caja.R, 0.0f);
    wall_poly[2].Set(-globalSetup->caja.R, globalSetup->caja.H); 
    wall_poly[3].Set(globalSetup->caja.R, globalSetup->caja.H);
    wall_poly[4].Set(globalSetup->caja.R, 2.0 * globalSetup->caja.r * sin(theta_rad));
    wall_poly[5].Set(globalSetup->caja.r * cos(theta_rad),
                     2.0 * globalSetup->caja.r * sin(theta_rad));

    b2ChainShape siloShape_poly;
    siloShape_poly.CreateChain(wall_poly, 6, wall_poly[0], wall_poly[5]);
    b2FixtureDef silo_fix;
    silo_fix.shape = &siloShape_poly;
    silo_fix.density = 0.0f;
    silo_fix.friction = globalSetup->caja.fric;
    silo->CreateFixture(&silo_fix);

    cout << "#\t- Silo creado." << endl;

    // Tapas	
    b2BodyDef tapa_piso;
    tapa_piso.position.Set(0.0f, 0.0f);
    tapa_piso.type = b2_staticBody;
    BodyData *tapaP = new BodyData;
    tapaP->isGrain = false;
    tapaP->nLados = 1;
    tapaP->gID = -110;
    tapa_piso.userData.pointer = uintptr_t(tapaP);
    b2Body *tapa_P = world.CreateBody(&tapa_piso);
    b2Vec2 v1(-globalSetup->caja.r * cos(theta_rad), 0.0f);
    b2Vec2 v2(globalSetup->caja.r * cos(theta_rad), 
              2.0 * globalSetup->caja.r * sin(theta_rad));
    b2EdgeShape tapa_p_f;
    tapa_p_f.SetTwoSided(v1, v2);
    b2FixtureDef tapa_p_Fix;
    tapa_p_Fix.shape = &tapa_p_f;
    tapa_p_Fix.density = 0.0f;
    tapa_p_Fix.friction = globalSetup->caja.fric;
    tapa_P->CreateFixture(&tapa_p_Fix);
    cout << "#\t- Tapa piso creada." << endl;
   
    // Generación de granos.
    float siloInf,siloSup,siloIzq,siloDer,x,y;
    float maxRadio = 0.0f;
    for (int i = 0; i < globalSetup->noTipoGranos; ++i) {
        if (globalSetup->granos[i]->radio > maxRadio) 
            maxRadio = globalSetup->granos[i]->radio;
    }
    float mg = 0.0;
    int noTotGranos = 0;
    siloInf = 2.7f * maxRadio + globalSetup->caja.y_max_hole;
    siloSup = globalSetup->caja.H - 1.5f * maxRadio;
    siloIzq = -globalSetup->caja.R + 1.7f * maxRadio;
    siloDer = globalSetup->caja.R - 1.7f * maxRadio;

    BodyData **gInfo;
    gInfo = new BodyData*[globalSetup->noTipoGranos];
    
    int contGid = 1;
    int *sumaTipo = new int[globalSetup->noTipoGranos] {};

    for (int i = 0; i < globalSetup->noTipoGranos; i++) { // Loop sobre tipos 
        gInfo[i] = new BodyData[globalSetup->granos[i]->noGranos]; // de granos
        for (int j = 0; j < globalSetup->granos[i]->noGranos; j++) { // Loop 
            // sobre el número de granos de cada tipo.
            x = rng.getAB(siloIzq, siloDer);
            y = rng.getAB(siloInf,siloSup);
            gInfo[i][j].tipo = i;
            gInfo[i][j].isGrain = true;
            gInfo[i][j].isIn = true;
            gInfo[i][j].nLados = globalSetup->granos[i]->nLados;
            gInfo[i][j].gID = contGid++;
            b2BodyDef bd;
            bd.type = b2_dynamicBody;
            bd.allowSleep = true;
            bd.bullet = globalSetup->isBullet;
            bd.position.Set(x,y);
            bd.angle = rng.getAB(-b2_pi,b2_pi);
            bd.userData.pointer = uintptr_t(&gInfo[i][j]);
            b2Body* grain = world.CreateBody(&bd);
            if (globalSetup->granos[i]->nLados == 1) {
                b2CircleShape circle;
                circle.m_radius = globalSetup->granos[i]->radio;
                b2FixtureDef fixDef;
                fixDef.shape = &circle;
                fixDef.density = globalSetup->granos[i]->dens;
                fixDef.friction = globalSetup->granos[i]->fric;
                fixDef.restitution = globalSetup->granos[i]->rest;
                grain->CreateFixture(&fixDef);
                mg += grain->GetMass();
                if (j == 0) cout << "#\t- Grano de tipo " << i 
                    << " creado con masa " << grain->GetMass() 
                        << " kg." << endl;
            }
            else {
                b2PolygonShape poly;
                int32 vertexCount = globalSetup->granos[i]->nLados;
                b2Vec2 vertices[8];
                for (int k = 0; k < globalSetup->granos[i]->nLados; k++) 
                    vertices[k].Set(globalSetup->granos[i]->vertices[k][0],
                            globalSetup->granos[i]->vertices[k][1]);
                poly.Set(vertices,vertexCount);
                b2FixtureDef fixDef;
                fixDef.shape = &poly;
                fixDef.density = globalSetup->granos[i]->dens;
                fixDef.friction = globalSetup->granos[i]->fric;
                fixDef.restitution = globalSetup->granos[i]->rest;
                grain->CreateFixture(&fixDef);
                mg += grain->GetMass();
                if (j == 0) cout << "#\t- Grano de tipo " << i + 1 
                    << " creado con masa " << grain->GetMass() 
                        << " kg." << endl;
            }
            noTotGranos++;
        } // Fin loop sobre el número de granos de cada tipo.
    } // Fin loop sobre tipo de granos	
    cout << "#\t- Número total de granos = " << noTotGranos << endl;
    cout << "#\t- Masa total de granos = " << mg << " kg."<< endl;
    float maxDistReiny = globalSetup->caja.r + 4.0f * maxRadio * 2.0f; // r + 4 diámetros de granos
    cout << "#\t- Límite de reinyección/eliminación = " << maxDistReiny << " m." << endl;
    float maxDistReiny_2 = maxDistReiny * maxDistReiny;
    float reinjectMinY = 0.0f;
    float reinjectMaxY = 0.0f;

    // Preparamos los parámetros de la simulación. 
    float timeStep = globalSetup->tStep;
    float timeS = 0.0;
    int32 pIterations = globalSetup->pIter;
    int32 vIterations = globalSetup->vIter;

    cout << "# Inicio de la simulación..." << endl;
    b2Vec2 pos; 
    float angle = 0.0, dist_out_2;
    int nGranosDesc = 0;
    BodyData *infGr;
    int deltaG = 0;


    // ---------------------------------------------------------
    // Variables para detección de avalanchas y atascos
    // ---------------------------------------------------------
    bool  inAvalanche = false;
    bool  inBlockage  = false;
    int   avalancheCount = 0;
    int   jamEventCount = 0;

    float lastParticleExitTime        = 0.0f;
    float avalancheStartTime          = 0.0f;
    float blockageStartTime           = 0.0f;
    int   exitedAtBlockageStart       = 0;    // contador al comenzar un atasco

    int   totalExitedParticles        = 0;   // total de granos descargados
    int   lastTotalExitedCount        = 0;   // valor en paso anterior
    int   particlesInCurrentAvalanche = 0;

    // Parámetros de tiempo (ajustables)
    const float BLOCKAGE_THRESHOLD      = 0.5f;   // s sin descargas → atasco
    const float STARTUP_GRACE          = 1.0f;   // no mirar atascos antes de este tiempo
    const float MIN_AVALANCHE_DURATION = 0.01f;  // descartar avalanchas ultra cortas

    // Archivo de salida para las avalanchas
    std::ofstream avalancheFile;
    std::string avalancheFileName = globalSetup->dirID + std::string("_avalanche_data.csv");
    avalancheFile.open(avalancheFileName.c_str());
    avalancheFile << "# index time_start time_end duration nparticles\n";
    avalancheFile.flush(); // escribir encabezado al instante



    // Llenado del silo y decantación
    bool saveFrm = (globalSetup->saveFrameFreq > 0 ? true : false);
    bool saveFlux = (globalSetup->fluxFreq > 0 ? true : false);
    int paso = 0;
    bool t_fill = (globalSetup->tFill > globalSetup->tStep ? true : false);
    // El tiempo de llenado funciona si los granos no son bullets
    cout << "# Llenado del silo y decantación... ";
    if (t_fill) cout << "(tiempo de espera " << globalSetup->tFill << " s.)" << endl;
    else cout << "(hasta el reposo)." << endl;
    do {
        world.Step(timeStep, pIterations, vIterations);
        timeS += timeStep;
        // Si es necesario, guardo el frame para graficar
         if (saveFrm && !(paso % globalSetup->saveFrameFreq)) {
             saveFrame(&world, paso, globalSetup);
         }
        paso++;
        if(t_fill && timeS > globalSetup->tFill) break;
        if(!t_fill && timeS > globalSetup->tMax) {
            cout << "Tiempo de llenado mayor que tiempo máximo de simulación. Fin." << endl;
            exit(1);
        }
    } while(isActive(&world));

    cout << "# ... tiempo de llenado (" << timeS << " s) finalizado." << endl;
    //saveFrame(&world, paso, globalSetup);

    // Altura máxima inicial de la columna granular para definir reinyección dinámica
    float initialTopY = 0.0f;
    for (b2Body *bd = world.GetBodyList(); bd; bd = bd->GetNext()) {
        BodyData* info = reinterpret_cast<BodyData*>(bd->GetUserData().pointer);
        if (!info || !info->isGrain || !info->isIn) continue;
        b2Vec2 p = bd->GetPosition();
        if (p.y > initialTopY) initialTopY = p.y;
    }
    reinjectMinY = initialTopY + 2.0f * static_cast<float>(globalSetup->Rmedio); // 2.5 diámetros por encima
    reinjectMaxY = reinjectMinY + 2.0f * static_cast<float>(globalSetup->Rmedio);
    //reinjectMinY = std::max(reinjectMinY, 0.75f * static_cast<float>(globalSetup->caja.H));
    //reinjectMaxY = std::max(reinjectMaxY, reinjectMinY + 2.0f * static_cast<float>(globalSetup->Rmedio));
    cout << "#\t- Altura inicial max: " << initialTopY << " m; reinyección en [" 
         << reinjectMinY << ", " << reinjectMaxY << "]" << endl;

    for (b2Body *bd = world.GetBodyList(); bd; bd = bd->GetNext()) {
        bd->SetAwake(true);
    }
                                
    //paso = 0;
    //saveFrame(&world, paso, globalSetup);
    float x_salida = 0.0f;
    float y_salida = 0.0f;
    world.DestroyBody(tapa_P);
    tapa_P = NULL;

    cout << "# Tapa removida, inicio de la descarga." << endl; 
    timeS = 0.0f;
    paso++;
    
    std::ofstream fileFlux;
    if (saveFlux) {
        fileFlux.open((globalSetup->fluxFile).c_str());
        fileFlux << "# theta = " << globalSetup->caja.theta << endl;
        fileFlux << "#grainDesc type time ";
        for (int i = 0; i < globalSetup->noTipoGranos; ++i) {
            fileFlux << "totalType_" << i + 1 << " ";
        }
        fileFlux << " Total" << endl;
    }
    // Guardo configuración inicial
    saveFrame(&world, paso, globalSetup);
    cout << "Frame " << paso << " guardado en " << timeS << endl;

    bool saveForces = (globalSetup->saveContactFreq > 0 ? true : false);
    bool saveCF;
    int forceFrame = 0; /*!< Contador de frames de fuerzas de contacto */

    // Inicio el avance del tiempo en la simulación
    for (;;) {
        if (globalSetup->maxAvalanchas > 0 &&
            avalancheCount >= globalSetup->maxAvalanchas) {
            cout << "# Número máximo de avalanchas alcanzado." << endl;
            break;
        }

        // Si es el momento guardo las fuerzas de contacto
        if (saveForces && timeS >= globalSetup->tStartCF && timeS <= globalSetup->tEndCF) {
            saveCF = true;
        }
        else saveCF = false;

        // Guardado de información de contacto
        if (saveCF && !(paso % globalSetup->saveContactFreq)) {
            forceFrame++;
            // Guardo info de las fuerzas de contacto
            saveContacts(&world, globalSetup->tStep, forceFrame, globalSetup);
            // Guardo info de los granos
            savePart(&world, forceFrame, globalSetup);
        }
        
        // Fuerza de rozamiento contra las paredes

        //b2Body* bd = world.GetBodyList();
        //while(bd) {
        //for (b2Body *bd = world.GetBodyList(); bd; bd = bd->GetNext()) {
            //infGr = (BodyData*)bd->GetUserData().pointer;
            //if(!infGr->isGrain) {
                ////bd = bd->GetNext();
                //continue;
            //}
            //if(!infGr->isIn) {
                ////bd = bd->GetNext();
                //continue;
            //}
            //if(bd->IsAwake()) {
               //vtmp = bd->GetLinearVelocity();
               //vtmpM = vtmp.Length();
               //if (vtmpM > vUmbral) {
                 //fric_force = globalSetup->caja.muBase * bd->GetMass();
                 //vtmp.Normalize();
                 //vtmp = -vtmp;
                 //vtmp *= fric_force;
                 //bd->ApplyForceToCenter(vtmp, true);
               //}
               //else {
                 //bd->SetLinearVelocity(b2Vec2(0.0f, 0.0f));     
               //}
            //}
            ////bd = bd->GetNext();
        //}

        // Cálculo y guardado de Q6
        if (globalSetup->q6_freq && !(paso % globalSetup->q6_freq)) {
            saveQ6(&world, paso, globalSetup);
        }

        // Cálculo y guardado del packing fraction
        //if (globalSetup->pf_freq && !(paso % globalSetup->pf_freq)) {
            //savepf(&world, paso, globalSetup);
        //}

        // Cálculo de descarga y reinyección
        deltaG = countDesc(&world, sumaTipo, paso, fileFlux, globalSetup);
        nGranosDesc += deltaG;

        // ---------------------------------------------------------
        // Lógica de flujo / avalanchas / atascos basada en countDesc
        // ---------------------------------------------------------
        {
            float currentTime = timeS;  // tiempo actual

            // Si hubo partículas descargadas en este paso
            if (deltaG > 0) {
                totalExitedParticles         += deltaG;
                lastParticleExitTime          = currentTime;


            }

            float timeSinceLastExit = currentTime - lastParticleExitTime;

            // Estado base: ni avalancha ni atasco
            if (!inAvalanche && !inBlockage) {

                // Si el total aumentó → comenzó avalancha
                if (totalExitedParticles > lastTotalExitedCount) {
                    inAvalanche = true;
                    avalancheStartTime          = currentTime;
                    particlesInCurrentAvalanche = deltaG;
                }
                // Si no aumentó y ya pasó el período de gracia → atasco inicial
                else if (currentTime > STARTUP_GRACE &&
                    timeSinceLastExit > BLOCKAGE_THRESHOLD) {
                    inBlockage        = true;
                    blockageStartTime = currentTime;
                    exitedAtBlockageStart = totalExitedParticles;
                    }
            }
            // Ya estamos en avalancha
            else if (inAvalanche) {

                // Acumulo descargas
                if (deltaG > 0) {
                    particlesInCurrentAvalanche += deltaG;
                }

                // Si pasa demasiado tiempo sin descargas → fin de avalancha y atasco
                if (timeSinceLastExit > BLOCKAGE_THRESHOLD) {
                    float avalDuration = currentTime - avalancheStartTime;

                    if (avalDuration >= MIN_AVALANCHE_DURATION &&
                        particlesInCurrentAvalanche > 0 &&
                        (globalSetup->maxAvalanchas <= 0 ||
                        avalancheCount < globalSetup->maxAvalanchas)) {

                        avalancheCount++;
                        avalancheFile << avalancheCount << " "
                                    << avalancheStartTime << " "
                                    << currentTime << " "
                                    << avalDuration << " "
                                    << particlesInCurrentAvalanche << "\n";
                        avalancheFile.flush();
                    }

                    inAvalanche = false;
                    inBlockage = true;
                    blockageStartTime = currentTime;
                    exitedAtBlockageStart = totalExitedParticles;
                    particlesInCurrentAvalanche = 0;
                }
            }
            // Estamos en atasco
            else if (inBlockage) {

                float blockageDuration = currentTime - blockageStartTime;

                // Si el atasco dura más que el tiempo configurado → intentar desatascar
                if (globalSetup->tiempoEsperaAtasco > 0.0 &&
                    blockageDuration >= globalSetup->tiempoEsperaAtasco) {

                // Registrar avalancha de tamaño cero si no hubo descargas
                if (totalExitedParticles == exitedAtBlockageStart &&
                    (globalSetup->maxAvalanchas <= 0 ||
                     avalancheCount < globalSetup->maxAvalanchas)) {

                    float avalDuration = currentTime - blockageStartTime;
                    if (avalDuration >= MIN_AVALANCHE_DURATION) {
                        avalancheCount++;
                        avalancheFile << avalancheCount << " "
                                      << blockageStartTime << " "
                                      << currentTime << " "
                                      << avalDuration << " "
                                      << 0 << "\n";
                        avalancheFile.flush();
                    }
                }

                std::cout << "# Atasco de " << blockageDuration
                 << " s, ejecutando unjam por raycast vertical." << std::endl;

                 // ---------------------------------------------------------
                // NUEVO: Guardar el frame del atasco antes de romperlo
                // ---------------------------------------------------------
                jamEventCount++; 
                saveAtascoFrame(&world, globalSetup, jamEventCount);
                // Contactos (entrada de presion_tensor_arch.py) y tensor de
                // estrés calculado en la simulación, del mismo instante
                saveAtascoContacts(&world, globalSetup, jamEventCount);
                saveAtascoTensor(&world, globalSetup, jamEventCount,
                                 currentTime);
                // ---------------------------------------------------------


                unjamArch(world, globalSetup, rng, siloIzq, siloDer, reinjectMinY, reinjectMaxY);

                // Reiniciamos el reloj de atasco después del intento
                blockageStartTime = currentTime;
                exitedAtBlockageStart = totalExitedParticles;
                    }

                    // Si vuelve a haber descarga → se reanuda el flujo
                    if (totalExitedParticles > lastTotalExitedCount) {
                        float blockageDur = currentTime - blockageStartTime;
                        std::cout << "# Flujo reanudado tras atasco de "
                        << blockageDur << " s" << std::endl;

                        inBlockage                  = false;
                        inAvalanche                 = true;
                        avalancheStartTime          = currentTime;
                        particlesInCurrentAvalanche = deltaG;
                    }
            }


            lastTotalExitedCount = totalExitedParticles;
        }




// --- BLOQUE MODIFICADO EN EL MAIN ---
        if (globalSetup->reinyect) {
            
            // 1. Preparamos la estructura (Lazy init).
            //    Aún no calculamos nada (computed = false).
            ReinjectionZones zones;
            zones.computed = false; 

            b2Body* bd = world.GetBodyList();
            while(bd){
                b2Body* b = bd;
                bd = bd->GetNext();
                
                infGr = (BodyData*)b->GetUserData().pointer;
                if(!infGr->isGrain) continue;
                if(infGr->isIn) continue;

                pos = b->GetPosition();
                dist_out_2 = (pos.x - x_salida) * (pos.x - x_salida) + 
                             (pos.y - y_salida) * (pos.y - y_salida);
                
                // Si no ha salido lo suficiente, saltamos
                if (dist_out_2 < maxDistReiny_2) continue;

                // Si llegamos aquí, hay que reinyectar o destruir
                if(globalSetup->reinyect) {
                    
                    // --- OPTIMIZACIÓN ---
                    // Si es la PRIMERA partícula que reinyectamos en este paso de tiempo,
                    // escaneamos el mundo para encontrar los huecos verticales.
                    if (!zones.computed) {
                        // Usamos maxRadio (que ya tienes definido arriba) para margenes seguros
                        float safeXMinMap = siloIzq + maxRadio; 
                        float safeXMaxMap = siloDer - maxRadio;
                        computeDensityMap(world, zones, safeXMinMap, safeXMaxMap, reinjectMinY, reinjectMaxY);
                    }
                    // --------------------

                    angle = b->GetAngle();
                    b2Fixture* f = b->GetFixtureList();
                    b2Shape* bs = (b2Shape*) f->GetShape();
                    float radio = bs->m_radius;

                    // 1. X Aleatoria (Distribuida horizontalmente)
                    float safeXMin = siloIzq + 2.0f*radio;
                    float safeXMax = siloDer - 2.0f*radio;
                    x = rng.getAB(safeXMin, safeXMax);

                    // 2. Y Inteligente (Basada en zonas de baja densidad pre-calculadas)
                    // Seleccionamos aleatoriamente uno de los bins vacíos encontrados
                    int randIdx = (int)rng.getAB(0, zones.validBins.size());
                    int selectedBin = zones.validBins[randIdx];
                    
                    float yBinMin = zones.minY + selectedBin * zones.binHeight;
                    float yBinMax = yBinMin + zones.binHeight;
                    
                    // Elegimos una altura aleatoria DENTRO de esa franja vacía
                    y = rng.getAB(yBinMin + radio*0.2f, yBinMax - radio*0.2f);

                    b2Vec2 newPos(x, y);
                    b->SetTransform(newPos, angle);

                    // Resetear física
                    b->SetLinearVelocity(b2Vec2(0.0f, 0.0f));
                    b->SetAngularVelocity(0.0f);
                    b->SetAwake(true);
                    infGr->isIn = true;
                }
                else {
                    world.DestroyBody(b);
                }
            }
        }

        world.Step(timeStep,pIterations,vIterations);
        paso++;
        timeS += timeStep;
        //cout << "paso: " << paso << " time: " << timeS << endl;

        // Si es necesario, guardo el frame para graficar
        if (saveFrm && !(paso % globalSetup->saveFrameFreq)) {
            saveFrame(&world, paso, globalSetup);
        }
    }
    // Fin del bucle principal; la salida ocurre solo al alcanzar maxAvalanchas.

    // ---------------------------------------------------------
    // Cierre: si la simulación termina en medio de una avalancha,
    // registramos esa avalancha parcial y cerramos el archivo.
    // ---------------------------------------------------------
    if (inAvalanche && particlesInCurrentAvalanche > 0) {
        float avalDuration = timeS - avalancheStartTime;

        if (avalDuration >= MIN_AVALANCHE_DURATION &&
            (globalSetup->maxAvalanchas <= 0 ||
            avalancheCount < globalSetup->maxAvalanchas)) {

            avalancheCount++;
            avalancheFile << avalancheCount << " "
                        << avalancheStartTime << " "
                        << timeS << " "
                        << avalDuration << " "
                        << particlesInCurrentAvalanche << "\n";
            avalancheFile.flush();
        }
    }

    if (avalancheFile.is_open()) {
        avalancheFile.close();
    }



    cout << "Descarga final: " << nGranosDesc << endl;
    for (int i = 0; i < globalSetup->noTipoGranos; ++i){
        cout << "Granos de tipo " << i + 1 << ": " << sumaTipo[i] << endl;
    }
    cout << "Simulación finalizada." << endl;

    
    // Clean up
    if (saveFlux) fileFlux.close();
    delete globalSetup;
    delete siloD;
    delete [] gInfo;
    delete [] sumaTipo;


    return 0;
}	
