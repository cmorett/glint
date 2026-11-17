// QA: deck.html?qa=1 — reports min font size, footer intrusions, off-slide elements, placeholder count
(function(){
  if(!new URLSearchParams(location.search).get('qa'))return;
  const out=[];let min=99,minEl='';
  document.querySelectorAll('.slide').forEach((sl,si)=>{
    const sr=sl.getBoundingClientRect();
    const walker=document.createTreeWalker(sl,NodeFilter.SHOW_TEXT);let n;
    while(n=walker.nextNode()){
      if(!n.textContent.trim())continue;const e=n.parentElement;
      let fs=parseFloat(getComputedStyle(e).fontSize);
      const svg=e.closest('svg');if(svg){const m=e.getScreenCTM?e.getScreenCTM():null;if(m)fs*=Math.hypot(m.a,m.b)}
      if(fs<min){min=fs;minEl=`s${si+1}: "${n.textContent.trim().slice(0,30)}"`}
      const r=document.createRange();r.selectNodeContents(n);const b=r.getBoundingClientRect();
      if(b.width===0)continue;
      const inFtr=e.closest('.ftr');
      if(!inFtr&&b.bottom-sr.top>1012)out.push(`s${si+1} FOOTER-INTRUSION "${n.textContent.trim().slice(0,40)}" bottom=${Math.round(b.bottom-sr.top)}`);
      if(b.right-sr.left>1920||b.left-sr.left<0||b.top-sr.top<0)out.push(`s${si+1} OFF-SLIDE "${n.textContent.trim().slice(0,40)}"`);
    }
    sl.querySelectorAll('*').forEach(e=>{if(e.closest('svg'))return;if(e.scrollWidth>e.clientWidth+1&&getComputedStyle(e).overflow==='hidden'&&!e.classList.contains('slide'))out.push(`s${si+1} CLIPPED <${e.tagName.toLowerCase()} class="${e.className}">`)});
  });
  out.unshift(`min font px: ${min.toFixed(1)} (${minEl})`,`placeholders: ${document.querySelectorAll('.ph').length}`,`slides: ${document.querySelectorAll('.slide').length}`);
  const pre=document.createElement('pre');pre.id='qa';pre.textContent='QA>>'+out.join('\nQA>>');document.body.appendChild(pre);
})();
