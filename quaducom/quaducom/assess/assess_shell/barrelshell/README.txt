Allgemeine Anmerkungen zum Export von 'geo_data' und 'state_data':
--> Die exportierten Dateien aus InfoCAD wurden �ber die Zwischenablage in einen Texteditor reinkkopiert. 
--> Die Namen stimmen mit den Nummern der Lastf�lle �berein.
--> Es wird Komma als Dezimaltrenner verwendet.

Zum Speichern der Werte in der Baumstrukturvon InfoCAD die folgenden Unterpunkte angeklicken:

### 'geo_data' ###

1) exportierte Datei: "Kontenkoordinaten.txt"
-->�ffnen unter:
"Strukturbeschreibung/Knotenkoordinaten" 
--> Header:
"Knotenkoordinaten [KNOTEN]
Nummer -- x[m] -- y[m] -- z[m]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

2) exportierte Datei: "Elementbeschreibung.txt"
--> �ffnen unter: 
"Strukturbeschreibung/Knotenkoordinaten" 
--> Header:
"Elementbeschreibung
Nummer -- Art -- Knoten-Nummer (1 .. 8) -- Querschn."
--> Anmerkungen:
hier bei Elementtyp SH46 nur 4 Knoten, Spalte 5-8 bleibt leer

3) exportierte Datei: "Querschnittswerte.txt"
-->�ffnen unter: 
"Ergebnisse/Protokolle/Finite Elemente" 
-->Werteblock unter der �berschrift "Querschnittswerte"
in eine Textdatei reinkopiert (nur Werte ohne �berschrift)
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

4) exportierte Datei: "Festhaltungen.txt"
"Strukturbeschreibung/Festhaltungen" 
--> Header:
"Festhaltungen
Knoten -- Drehung des Lagersystem um Achse (x,y,z) [�] -- ux,y,z -- phix,y,z -- Zugausfall"
--> Anmerkungen:
Die Drehwinkel geben an, um wie viel Grad das Lager-KO-Sytsem gegen�ber den globalen Achsen gedreht wurde.


### 'state_data' ###

Zum Export der Schnittgr��en in den Elementmittelpunkten unter dem Men�punkt "Berechnung/Einstellungen"
die Option "Berechnungsort Fl�chenschnittgr��en" auf "Schwerpunkt" einstellen und Berechnung starten.

4) exportierte Dateien: "LC1.txt", "LC2.txt", "LC3.txt", (...) im Ordner "Flaechenschnittgroessen/Schwerpunkt" 
-->�ffnen unter: 
"Ergebnisse/Schnittgr��en Fl�chenelemente" 
-->Header:
"Schnittgr��en Fl�chenelemente:
Element -- nx [kN/m] -- ny [kN/m] -- nxy [kN/m] -- mx [kNm/m] -- my [kNm/m] -- mxy [kNm/m]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

5) exportierte Dateien: "LC1.txt", "LC2.txt", "LC3.txt", (...) im Ordner "Auflagerreaktionen.txt"
-->�ffnen unter: 
"Ergebnisse/Auflagerreaktionen" 
-->Header:
"Auflagerreaktionen [AUFLR.1]:
Knoten -- Rx [kN] -- Ry [kN] -- Rz [kN] -- Mx [kNm] -- My [kNm] -- Mz [kNm]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

6) exportierte Dateien: "LC1.txt", "LC2.txt", "LC3.txt", (...) im Ordner "Knotendeformationen.txt"
-->�ffnen unter: 
"Ergebnisse/Knotendeformationen" 
-->Header:
"Knotendeformationen LF# [DEFORM.#]:
Knoten -- ux [m] -- uy [m -- uz [m] -- phi.x [rad] -- phi.y [rad] -- phi.z [rad]"
-->mit Hilfe des Texteditors in der Datei "," durch "." ersetzen

