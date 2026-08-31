# Metodologia tecnico-cientifica del modelo de derrame superficial de hidrocarburos

## 1. Objetivo y alcance

El modelo simula la evolucion espacio-temporal de un derrame superficial sobre terreno, con enfasis en hidrodinamica de lamina somera, interaccion con topografia e intercambio vertical con el suelo. El enfoque permite cuantificar:

1. Huella de afectacion y avance del frente de derrame.
2. Volumen superficial remanente en el tiempo.
3. Perdidas por infiltracion y por atenuacion ambiental simplificada.
4. Aportes externos por lluvia y por descarga puntual.
5. Cierre de balance de masa para control de calidad tecnico.

## 2. Marco conceptual

La metodologia representa el hidrocarburo como una lamina de espesor pequeno respecto de las escalas horizontales del terreno. Bajo esta hipotesis, el movimiento se describe mediante ecuaciones de aguas someras bidimensionales, con forzamiento gravitacional, efectos de pendiente topografica y disipacion por friccion basal.

Se incorpora ademas un esquema de infiltracion tipo Green-Ampt y terminos de sumidero para evaporacion y degradacion, de forma que el modelo cierre el balance volumetrico de manera trazable.

## 3. Formulacion fisico-matematica

Variables de estado:

- $h(x,y,t)$: espesor de la lamina superficial.
- $q_x = h u$: descarga unitaria en direccion $x$.
- $q_y = h v$: descarga unitaria en direccion $y$.
- $z(x,y)$: elevacion del terreno.
- $g$: aceleracion de gravedad.

Ecuacion de continuidad:

$$
\frac{\partial h}{\partial t} + \frac{\partial q_x}{\partial x} + \frac{\partial q_y}{\partial y} = R + Q_p - I - E
$$

Ecuaciones de momento:

$$
\frac{\partial q_x}{\partial t} + \frac{\partial}{\partial x}\left(\frac{q_x^2}{h} + \frac{1}{2}gh^2\right) + \frac{\partial}{\partial y}\left(\frac{q_x q_y}{h}\right) = -gh\frac{\partial z}{\partial x} - S_{fx}
$$

$$
\frac{\partial q_y}{\partial t} + \frac{\partial}{\partial x}\left(\frac{q_x q_y}{h}\right) + \frac{\partial}{\partial y}\left(\frac{q_y^2}{h} + \frac{1}{2}gh^2\right) = -gh\frac{\partial z}{\partial y} - S_{fy}
$$

Donde:

- $R$ representa lluvia uniforme.
- $Q_p$ representa inyeccion puntual externa.
- $I$ representa infiltracion al suelo.
- $E$ representa perdida por evaporacion y degradacion.

## 4. Parametrizacion de procesos

### 4.1 Friccion hidraulica

Se consideran dos cierres alternativos de friccion basal:

1. Manning, apropiado para rugosidad equivalente de superficie.
2. Darcy-Weisbach, apropiado para representar perdida de energia mediante factor de friccion.

Ambos cierres disipan momento y controlan la velocidad de propagacion del frente.

### 4.2 Infiltracion

La infiltracion se representa con Green-Ampt:

$$
i = K_s\left(1 + \frac{\psi_f \Delta\theta}{F}\right)
$$

con:

- $K_s$: conductividad hidraulica saturada.
- $\psi_f$: succion capilar efectiva.
- $\Delta\theta$: deficit de porosidad.
- $F$: infiltracion acumulada.

Se aplica regularizacion en $F$ para evitar tasas iniciales singulares y se limita la infiltracion al espesor disponible en superficie.

### 4.3 Atenuacion superficial

Evaporacion y degradacion se modelan como un sumidero agregado con tasa areal constante. Este termino es util para analisis de escenario y priorizacion de respuesta.

### 4.4 Forzantes externos

- Lluvia: aporte uniforme de espesor en el dominio.
- Fuente puntual: aporte localizado dependiente del tiempo.

## 5. Estrategia numerica

La resolucion emplea un esquema de volumen finito en malla cartesiana regular, con las siguientes caracteristicas:

1. Flujo numerico de Rusanov para resolver interfaces y capturar frentes de propagacion de forma robusta.
2. Reconstruccion hidrostatica para tratamiento consistente de topografia y preservacion del equilibrio de reposo.
3. Correccion de positividad para impedir espesores negativos.
4. Tratamiento de celdas secas para estabilidad en zonas de frente humedo-seco.
5. Condiciones de borde reflectivas (pared impermeable), adecuadas para dominios cerrados.

## 6. Estabilidad temporal

El paso de tiempo es adaptativo con criterio CFL:

$$
\Delta t = \mathrm{CFL}\,\frac{\min(\Delta x,\Delta y)}{\max\left(|u| + \sqrt{gh},\; |v| + \sqrt{gh}\right)}
$$

Este control garantiza estabilidad numerica y coherencia fisica en la propagacion del derrame.

## 7. Control de balance de masa

El marco metodologico calcula y reporta:

1. Masa/volumen superficial inicial.
2. Entradas acumuladas por lluvia.
3. Entradas acumuladas por fuente puntual.
4. Masa/volumen infiltrado acumulado.
5. Masa/volumen perdido por evaporacion y degradacion.
6. Masa/volumen superficial final.
7. Error de cierre porcentual del balance.

Este esquema de auditoria es clave para trazabilidad frente a operador y autoridad ambiental.

## 8. Aplicacion a evaluacion de riesgo ambiental

La metodologia permite derivar productos utiles para gestion de riesgo:

1. Mapas de espesor del hidrocarburo por tiempo.
2. Delimitacion de huella de afectacion y direccion preferente de propagacion.
3. Estimacion de persistencia superficial y remocion natural.
4. Priorizacion de zonas para contencion, recuperacion y remediacion.

## 9. Supuestos, alcances y limites

1. El hidrocarburo se trata como lamina equivalente monofasica (enfoque volumetrico superficial).
2. No se representa en detalle la composicion multicomponente ni la termodinamica del producto.
3. No se explicita transporte en zona no saturada profunda.
4. Evaporacion y degradacion se simplifican como tasa agregada.
5. La precision depende de calidad de topografia, parametrizacion de friccion e infiltracion, y condicion inicial del derrame.

## 10. Uso recomendado en ingenieria

El enfoque es adecuado para prefactibilidad, analisis de escenarios y soporte a decisiones de respuesta inicial. Para aplicaciones regulatorias o de diseno de detalle se recomienda complementar con:

1. Calibracion con datos de campo o eventos historicos.
2. Analisis de sensibilidad e incertidumbre de parametros.
3. Escenarios estacionales de lluvia y humedad antecedente.
4. Conversión volumen-masa con propiedades fisicas especificas del hidrocarburo evaluado.
