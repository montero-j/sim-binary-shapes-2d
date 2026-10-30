/*! \file jam_output.hpp
 * \brief Salidas de cada atasco: contactos y tensor de estrés por grano.
 *
 * Se llaman junto con saveAtascoFrame, antes de romper el arco. En ese
 * momento las posiciones y los impulsos de los manifolds corresponden al
 * mismo Step (el último), así que las tres salidas son consistentes.
 *
 * Convención de Box2D 2.4 (b2ContactSolver): la normal n del manifold
 * apunta del cuerpo A al cuerpo B y la tangente es t = b2Cross(n, 1) =
 * (n.y, -n.x). La fuerza media durante el paso sobre B, en el punto i, es
 *
 *     F_B = (normalImpulse n + tangentImpulse t) / dt,    F_A = -F_B.
 *
 * Tensor de estrés de cada grano p (fórmula de Love, compresión < 0):
 *
 *     s_ij = (1/A_p) sum_c F_i l_j,   l = punto de contacto - centro de masa,
 *
 * con A_p el área del grano según Box2D (la que da su masa: para polígonos,
 * el área de los vértices, sin la piel b2_polygonRadius). sn es la parte
 * debida solo a las fuerzas normales; la tangencial es s - sn. Ver
 * docs/stress/poligonos.tex en unav-walker para las propiedades con
 * polígonos (sn no simétrico, la fricción aporta a la traza).
 */

#ifndef _JAM_OUTPUT_H
#define _JAM_OUTPUT_H

#include <box2d/box2d.h>
#include "globalSetup.hpp"

/*! Escribe frames_<dirID>/atasco_contacts_<jam>.dat con el formato del
 * binario de producción: una fila por punto de manifold de cada contacto
 * que se toca (granos y paredes), "ID_A ID_B PosX PosY NormalImpulse
 * TangentImpulse", con impulsos (no fuerzas) y la normal de A hacia B.
 * Es la entrada de analysis/3_presion_y_mapas/presion_tensor_arch.py. */
void saveAtascoContacts(b2World *world, const GlobalSetup *gs, int jamIndex);

/*! Calcula en la simulación el tensor de estrés de cada grano y escribe:
 *  - frames_<dirID>/atasco_tensor_<jam>.dat: una fila por grano activo
 *    (isGrain e isIn, los mismos del atasco_<jam>.xy);
 *  - frames_<dirID>/atasco_contacts_full_<jam>.dat: una fila por punto de
 *    manifold, con la normal de Box2D, las fuerzas y los centros, para
 *    comparar contacto por contacto con la normal que infiere el análisis.
 */
void saveAtascoTensor(b2World *world, const GlobalSetup *gs, int jamIndex,
                      double tSim);

#endif
