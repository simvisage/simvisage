__author__ = 'campsb'
__date__ = 'Nov 8, 2017'


import os

import matplotlib.pyplot as plt     # Bibliothek fuer das erstellen von Diagrammen
from matresdev.db.simdb import SimDB
import numpy as np                  # Bibliothek fuer grosse Datenmengen
import pandas as pd                 # Bibliothek fuer grosse tabellarische Daten
simdb = SimDB()


# define the path to the file
test_file_path = os.path.join(simdb.exdata_dir,
                              'compression_tests', 'cylinder_tests',
                              )

test_file = 'WinConFat_CT_120_1_5.csv'
rohdaten_dir = os.path.join(test_file_path, test_file)


##########################################################################
# Variablen
# Anzahl der Zeilen fuer die die Berechnungen durchgefuehrt werden soll
nrows = None
# fuer alle gilt nrows=None
versuch = 'CT_120_1_15'

script_dir = os.path.dirname(__file__)
# Dateipfad des Scriptes
#rohdaten_dir = os.path.join(script_dir, 'rohdaten_CT', '%s.csv' % (versuch))
# Dateipfad der Rohdaten
diagramm_dir = os.path.join(
    script_dir, 'auswertung_CT', 'diagramme_%s' % (versuch))
# Speicherort der Diagramme

if not os.path.isdir(diagramm_dir):
    # Erstellt einen Ordner mit Name=Versuch, falls noch nicht vorhanden
    os.makedirs(diagramm_dir)


##########################################################################
# Erscheinungsbild Diagramme

font = {'name': 'arial', 'size': 8}    # Schriftart und Schriftgroesse

line_widths = {'main': 1.5, 'sec': 0.5, 'graph': 1.0}
# Linienstaerken Haupt-(main), Nebenachsen (sec) und Kurve (Graph)

ticks = {'width_major': 1.5, 'width_minor': 0.5,
         'size_major': 3.5, 'size_minor': 3.0, 'log': False}
# Linienstaerke der Achsenbemassung (Ticks) und logarithmische Darstellung

font_title = {'size': 10, 'weight': 'bold'}
# Schriftgroesse und -breite

font_axistitle = {'size': 8, 'weight': 'bold'}
# Schriftgroesse und -breite Achsenbeschriftung

fig_dimension = {'length': 6.30, 'hight': 4.72, 'dpi': 300}
# Breite, Hoehe und Aufloesung der Diagramme

adjust_graph = {'bottom': 0.21, 'top': 0.9, 'right': 0.95, 'left': 0.18}
# Position der Kurve

adjust_legend = {'size': 0.0005, 'spacing': 0.2}
# Schriftgroesse und Abstand Legende
# size=0.001: Legende nicht sichtbar, sonst 6

adjust_xlim = {'xWA_min': 0., 'xWA_max': 1.5}
# Wertebereich X-Achse
adjust_ylim = {'yF_min': 0., 'yF_max': 1000}
# Wertebereich Y-Achse


##########################################################################
# Funktionsdefinition

plt.rc("font", family=font['name'], size=font['size'])
# Definition der Schriftgroesse und -name


def konvertiere_kraft_positiv(kraft):
    return kraft * -1.0               # Kraft wird umgerechnet


def definiere_design(ax, title, xlabel, ylabel, xlims, ylims, log=ticks['log']):
                                    # Funktion fuer Design der Diagramme
    # def mittelwert_WA(WA_1, WA_2, WA_3):
    # return (WA_1+WA_2+WA_3)/3

    ax.set_title(title, size=font_title['size'],
                 fontweight=font_title['weight'])
    # Diagrammtitel
    if log:
        ax.set_xscale('log')        # logarithmische Darstellung

    ax.set_xlim(xlims[0], xlims[1])  # Achsenabschnitte
    ax.set_ylim(ylims[0], ylims[1])

    ax.set_xlabel(
        xlabel, size=font_axistitle['size'], fontweight=font_axistitle['weight'])
    # X-Achsenbeschriftung
    ax.set_ylabel(
        ylabel, size=font_axistitle['size'], fontweight=font_axistitle['weight'])
    # Y-Achsenbeschriftung
    plt.subplots_adjust(bottom=adjust_graph['bottom'], top=adjust_graph['top'],
                        left=adjust_graph['left'], right=adjust_graph['right'])
    # Ausrichten des Plots
    ax.tick_params(axis='both', which='major', direction='out',
                   size=ticks['size_major'], width=ticks['width_major'])
    ax.tick_params(axis='both', which='minor', direction='out',
                   size=ticks['size_minor'], width=ticks['width_minor'])
    # Ticks entfernen
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()      # nicht gebrauchte Ticks werden entfernt

    ax.spines["left"].set_linewidth(line_widths['main'])
    ax.spines["bottom"].set_linewidth(line_widths['main'])
    ax.spines["right"].set_linewidth(line_widths['sec'])
    ax.spines["top"].set_linewidth(line_widths['sec'])
    ax.xaxis.get_major_formatter().set_powerlimits((0, 0))
    # exponentielle Darstellung


##########################################################################
# Erstellung der Datenbasis fuer Diagramme

names = ['Zeit [s]', 'Kraft F [kN]', 'Maschinenweg [mm]',
         'WA_1 [mm]', 'WA_2 [mm]', 'WA_3 [mm]']
# Definition der Spaltennamen
F = names[1]
#WA_i='MW WA'
WA_1 = names[3]
WA_2 = names[4]
# Definition der Abkuerzungen alias Spalten
WA_3 = names[5]

daten = pd.read_csv(rohdaten_dir, sep=';', decimal=',',
                    encoding='iso-8859-1', skiprows=2, names=names, nrows=nrows)
# ...read_csv(Pfad, Spaltentrennung, Dezimalzeichen, Codierung, uebersprungene Zeilen, Spaltennamen, Anzahl Zeilen)


##########################################################################
# Anwendung der Funktionen

daten[F] = daten[F].map(konvertiere_kraft_positiv)
# .map fuehrt Funktion in jeder Zeile der Spalte [F] aus und ueberschreibt daten[F] mit positiven Zahlenwerten


##########################################################################
# Erstellung der Daten zum Plotten und Loeschung der leeren Zeilen

daten_cut = daten[0:daten[F].idxmax() + 1]
# alle Daten nach Fmax werden verworfen

WA1_daten = daten_cut.ix[:, [F, WA_1]]     # Erstellt Datensatz
# daten.ix erstellt einen Index
# : definiert alle zeilen, F die Spalten
WA1_daten.dropna(axis=0, inplace=True)
# loescht alle leeren Zeilen, in denen kein Weg vorhanden ist

WA2_daten = daten_cut.ix[:, [F, WA_2]]     # Erstellt Datensatz
# daten.ix erstellt einen Index
# : definiert alle zeilen, F die Spalten
WA2_daten.dropna(axis=0, inplace=True)
# loescht alle leeren Zeilen, in denen kein Weg vorhanden ist

WA3_daten = daten_cut.ix[:, [F, WA_3]]     # Erstellt Datensatz
# daten.ix erstellt einen Index
# : definiert alle zeilen, F die Spalten
WA3_daten.dropna(axis=0, inplace=True)
# loescht alle leeren Zeilen, in denen kein Weg vorhanden ist

##########################################################################
# Plot WA_1

fig = plt.figure(figsize=(
    fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])
# Figure zum Plotten wird im Vorfeld angelegt

ax1 = fig.add_subplot(111)

WA1_daten.plot(x=WA_1, y=F, ax=ax1, color=(
    '0.0'), linewidth=line_widths['graph'])
# datenname.plot(X-Achse, Y-Achse, Achsendef., Farbe (1.0=Weiss,
# 0.0=Schwarz, Linienstaerke)
definiere_design(ax1, versuch, 'Verformung w [mm]', 'Kraft F [kN]', (
    adjust_xlim['xWA_min'], adjust_xlim['xWA_max']), (adjust_ylim['yF_min'], adjust_ylim['yF_max']))
# Achsenbeschriftung (nach.plot-Funktion)
ax1.legend(loc=2, prop={
           'size': adjust_legend['size']}, labelspacing=adjust_legend['spacing'])
# Position Legende (1=oben rechts, 2=obenlinks, 3=unten links, 4=unten rechts)
fig.savefig(diagramm_dir + '//WA_1.eps')


# Plot WA_2

fig = plt.figure(figsize=(
    fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])
# Figure zum Plotten wird im Vorfeld angelegt

ax2 = fig.add_subplot(111)

WA2_daten.plot(x=WA_2, y=F, ax=ax2, color=(
    '0.0'), linewidth=line_widths['graph'])
# datenname.plot(X-Achse, Y-Achse, Achsendef., Farbe (1.0=Weiss,
# 0.0=Schwarz, Linienstaerke)
definiere_design(ax2, versuch, 'Verformung w [mm]', 'Kraft F [kN]', (
    adjust_xlim['xWA_min'], adjust_xlim['xWA_max']), (adjust_ylim['yF_min'], adjust_ylim['yF_max']))
# Achsenbeschriftung (nach.plot-Funktion)
ax2.legend(loc=2, prop={
           'size': adjust_legend['size']}, labelspacing=adjust_legend['spacing'])
# Position Legende (1=oben rechts, 2=obenlinks, 3=unten links, 4=unten rechts)
fig.savefig(diagramm_dir + '//WA_2.eps')


# Plot WA_3

fig = plt.figure(figsize=(
    fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])
# Figure zum Plotten wird im Vorfeld angelegt

ax3 = fig.add_subplot(111)

WA3_daten.plot(x=WA_3, y=F, ax=ax3, color=(
    '0.0'), linewidth=line_widths['graph'])
# datenname.plot(X-Achse, Y-Achse, Achsendef., Farbe (1.0=Weiss,
# 0.0=Schwarz, Linienstaerke)
definiere_design(ax3, versuch, 'Verformung w [mm]', 'Kraft F [kN]', (
    adjust_xlim['xWA_min'], adjust_xlim['xWA_max']), (adjust_ylim['yF_min'], adjust_ylim['yF_max']))
# Achsenbeschriftung (nach.plot-Funktion)
ax3.legend(loc=2, prop={
           'size': adjust_legend['size']}, labelspacing=adjust_legend['spacing'])
# Position Legende (1=oben rechts, 2=obenlinks, 3=unten links, 4=unten rechts)
fig.savefig(diagramm_dir + '//WA_3.eps')

# Plot alle WA_i

fig = plt.figure(figsize=(
    fig_dimension['length'], fig_dimension['hight']), dpi=fig_dimension['dpi'])
# Figure zum Plotten wird im Vorfeld angelegt

ax1 = fig.add_subplot(111)

WA1_daten.plot(x=WA_1, y=F, ax=ax1, color=(
    '0.0'), linewidth=line_widths['graph'])
WA2_daten.plot(x=WA_2, y=F, ax=ax1, color=(
    '0.4'), linewidth=line_widths['graph'])
WA3_daten.plot(x=WA_3, y=F, ax=ax1, color=(
    '0.8'), linewidth=line_widths['graph'])
# datenname.plot(X-Achse, Y-Achse, Achsendef., Farbe (1.0=Weiss,
# 0.0=Schwarz, Linienstaerke)
definiere_design(ax1, versuch, 'Verformung w [mm]', 'Kraft F [kN]', (
    adjust_xlim['xWA_min'], adjust_xlim['xWA_max']), (adjust_ylim['yF_min'], adjust_ylim['yF_max']))
# Achsenbeschriftung (nach.plot-Funktion)
ax1.legend(loc=2, prop={
           'size': adjust_legend['size']}, labelspacing=adjust_legend['spacing'])
# Position Legende (1=oben rechts, 2=obenlinks, 3=unten links, 4=unten rechts)
plt.plot()

plt.show()
#fig.savefig(diagramm_dir + '//WA_alle.eps')

##########################################################################
# Ausgabe der maximalen Werte fuer statische Versuche

print((daten.iloc[daten[F].idxmax()]))
