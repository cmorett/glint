#!/usr/bin/env bash
# Renders deck.html -> slide-N.png + Glint_YPrize_Phase1.pdf, and brand.html -> brand.png
# Uses the installed Chrome in headless mode (Playwright fallback).
# Usage: ./render.sh [all|slides|pdf|brand|N]
set -e
DIR="$(cd "$(dirname "$0")" && pwd)"
C="${CHROME:-/c/Program Files/Google/Chrome/Application/chrome.exe}"
WD="$(cygpath -m "$DIR")"; URL="file:///$(cygpath -m "$DIR")"
FLAGS=(--headless=new --disable-gpu --hide-scrollbars --force-device-scale-factor=1 --virtual-time-budget=4000)
shot(){ "$C" "${FLAGS[@]}" --window-size=$3 --screenshot="$WD/$2" "$URL/$1" 2>&1 | grep -v external_registry || true; }
what="${1:-all}"
if [[ "$what" == all || "$what" == brand ]]; then shot "brand.html" brand.png 1920,3640; fi
if [[ "$what" == all || "$what" == slides ]]; then for i in 1 2 3 4 5; do shot "deck.html?s=$i" slide-$i.png 1920,1080; done; fi
if [[ "$what" =~ ^[1-5]$ ]]; then shot "deck.html?s=$what" slide-$what.png 1920,1080; fi
if [[ "$what" == all || "$what" == pdf ]]; then
  "$C" "${FLAGS[@]}" --no-pdf-header-footer --print-to-pdf="$WD/Glint_YPrize_Phase1.pdf" "$URL/deck.html" >/dev/null 2>&1; echo "wrote Glint_YPrize_Phase1.pdf"
fi
