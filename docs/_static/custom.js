document.addEventListener("DOMContentLoaded", function () {
    document.querySelectorAll(".collapsible-toc li").forEach(function (item) {
        const sublist = item.querySelector("ul");
        const link = item.querySelector("a");
        if (sublist && link) {
            const toggle = document.createElement("span");
            toggle.innerHTML = `
                <svg class="toctree-expander" viewBox="0 0 24 24" width="16" height="16" stroke="currentColor" fill="none" stroke-width="2">
                    <polyline points="6 9 12 15 18 9"></polyline>
                </svg>
            `;
            toggle.style.cursor = "pointer";
            toggle.style.marginLeft = "8px";

            toggle.addEventListener("click", function () {
                if (sublist.style.display === "none") {
                    sublist.style.display = "block";
                    toggle.querySelector("polyline").setAttribute("points", "6 15 12 9 18 15");
                } else {
                    sublist.style.display = "none";
                    toggle.querySelector("polyline").setAttribute("points", "6 9 12 15 18 9");
                }
            });

            sublist.style.display = "none";
            link.insertAdjacentElement("afterend", toggle);
        }
    });
});
