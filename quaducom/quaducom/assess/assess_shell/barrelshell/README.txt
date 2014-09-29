Allgemeine Anmerkungen zum Export von 'geo_data' und 'state_data':
--> Die exportierten Dateien aus InfoCAD wurden über die Zwischenablage in einen Texteditor reinkkopiert. 
--> Die Namen stimmen mit den Nummern der Lastfälle überein.
--> Es wird Komma als Dezimaltrenner verwendet.

Zum Speichern der Werte in der Baumstrukturvon InfoCAD die folgenden Unterpunkte angeklicken:

### 'geo_data' ###

1) exportierte Datei: "Kontenkoordinaten.txt"
-->Öffnen unter:
"Strukturbeschreibung/Knotenkoordinaten" 
--> Header:
"Knotenkoordinaten [KNOTEN]
Nummer -- x[m] -- y[m] -- z[m]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

2) exportierte Datei: "Elementbeschreibung.txt"
--> Öffnen unter: 
"Strukturbeschreibung/Knotenkoordinaten" 
--> Header:
"Elementbeschreibung
Nummer -- Art -- Knoten-Nummer (1 .. 8) -- Querschn."
--> Anmerkungen:
hier bei Elementtyp SH46 nur 4 Knoten, Spalte 5-8 bleibt leer

3) exportierte Datei: "Querschnittswerte.txt"
-->Öffnen unter: 
"Ergebnisse/Protokolle/Finite Elemente" 
-->Werteblock unter der Überschrift "Querschnittswerte"
in eine Textdatei reinkopiert (nur Werte ohne Überschrift)
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

4) exportierte Datei: "Festhaltungen.txt"
"Strukturbeschreibung/Festhaltungen" 
--> Header:
"Festhaltungen
Knoten -- Drehung des Lagersystem um Achse (x,y,z) [°] -- ux,y,z -- phix,y,z -- Zugausfall"
--> Anmerkungen:
Die Drehwinkel geben an, um wie viel Grad das Lager-KO-Sytsem gegenüber den globalen Achsen gedreht wurde.


### 'state_data' ###

Zum Export der Schnittgrößen in den Elementmittelpunkten unter dem Menüpunkt "Berechnung/Einstellungen"
die Option "Berechnungsort Flächenschnittgrößen" auf "Schwerpunkt" einstellen und Berechnung starten.

4) exportierte Dateien: "LC1.txt", "LC2.txt", "LC3.txt", (...) im Ordner "Flaechenschnittgroessen/Schwerpunkt" 
-->Öffnen unter: 
"Ergebnisse/Schnittgrößen Flächenelemente" 
-->Header:
"Schnittgrößen Flächenelemente:
Element -- nx [kN/m] -- ny [kN/m] -- nxy [kN/m] -- mx [kNm/m] -- my [kNm/m] -- mxy [kNm/m]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

5) exportierte Dateien: "LC1.txt", "LC2.txt", "LC3.txt", (...) im Ordner "Auflagerreaktionen.txt"
-->Öffnen unter: 
"Ergebnisse/Auflagerreaktionen" 
-->Header:
"Auflagerreaktionen [AUFLR.1]:
Knoten -- Rx [kN] -- Ry [kN] -- Rz [kN] -- Mx [kNm] -- My [kNm] -- Mz [kNm]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

6) exportierte Dateien: "LC1.txt", "LC2.txt", "LC3.txt", (...) im Ordner "Knotendeformationen.txt"
-->Öffnen unter: 
"Ergebnisse/Knotendeformationen" 
-->Header:
"Knotendeformationen LF# [DEFORM.#]:
Knoten -- ux [m] -- uy [m -- uz [m] -- phi.x [rad] -- phi.y [rad] -- phi.z [rad]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

