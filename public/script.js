
// Fonction pour envoyer la requête Prolog au serveur via une requête HTTP POST
async function sendQuery() {
    const queryText = 'exexe.';

    try {
        const response = await fetch('/execute-prolog', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({ query: queryText })
        });

        if (!response.ok) {
            throw new Error('Erreur lors de la requête Prolog');




        }

    // Récupérer les résultats de la requête Prolog
        const data = await response.json();
        const results = data.results;

        // Afficher les résultats dans la zone de sortie
       // displayResults(results);



        // Récupérer le contenu du fichier output.txt et l'afficher dans la zone de sortie
        fetchOutput();
    } catch (error) {
        console.error('Erreur lors de la requête Prolog:', error);
        const outputDiv = document.getElementById('output');
        outputDiv.innerHTML = ''; // Effacer le contenu précédent
        const errorMessage = document.createElement('p');
        errorMessage.textContent = 'Erreur lors de l\'exécution de la requête Prolog.';
        outputDiv.appendChild(errorMessage);
    }
}
// Afficher le script 
function afficherScript() {
    // Envoyer une requête au serveur pour récupérer le contenu des fichiers
    fetch('/get-files-content')
    .then(response => response.json())
    .then(data => {
        // Ouvrir une nouvelle fenêtre et afficher le contenu des fichiers
        const newWindow = window.open('', '_blank');
        newWindow.document.write('<pre>Contenu de input.txt :</pre>');
        newWindow.document.write(`<pre>${data.inputContent}</pre>`);
        newWindow.document.write('<pre>Contenu de getg.txt :</pre>');
        newWindow.document.write(`<pre>${data.getgContent}</pre>`);
    })
    .catch(error => console.error('Erreur lors de la récupération du contenu des fichiers :', error));
}




/*// Fonction pour afficher les résultats dans la zone de sortie de l'interface web
function displayResults(results) {
    const outputDiv = document.getElementById('output');
    outputDiv.innerHTML = ''; // Effacer le contenu précédent

    if (results && results.length > 0) {
        // Afficher chaque résultat dans un élément <p>
        results.forEach((result, index) => {
            const resultText = `Solution ${index + 1}: ${result}`;
            const resultElement = document.createElement('p');
            resultElement.textContent = resultText;
            outputDiv.appendChild(resultElement);
        });
    } else {
        // Aucun résultat trouvé ou résultats non valides
        const errorMessage = document.createElement('p');
        errorMessage.textContent = 'Aucun résultat valide trouvé.';
        outputDiv.appendChild(errorMessage);
    }
}*/



// Fonction pour récupérer le contenu du fichier output.txt et l'afficher dans la zone de sortie
async function fetchOutput() {
    try {
        const response = await fetch('/get-output');
        if (!response.ok) {
            throw new Error('Erreur lors de la récupération du contenu du fichier output.txt');
        }
        const output = await response.text();
        const outputDiv = document.getElementById('output');
        outputDiv.innerHTML += `<p></p><pre>${output}</pre>`;
    } catch (error) {
        console.error('Erreur lors de la récupération du contenu du fichier output.txt:', error);
    }
}


function afficherCodeExemples() {
    // Rediriger vers ecran.html
    window.location.href = 'code.html';
}
function afficherScript(){
    window.location.href ='script.html';
}





// Fonction générique pour renommer un fichier d'entrée spécifique
function renameFileToInput(endpoint) {
    fetch(endpoint, {
        method: 'POST'
    })
    .then(response => response.text())
    .then(data => {
        console.log(data);
        // Affiche un message de succès ou effectue d'autres actions si nécessaire
    })
    .catch(error => {
        console.error('Erreur lors du renommage du fichier:', error);
    });
}

// Fonctions pour chaque bouton spécifique (but1, but2, but3, but4, but5)
function renameBut1ToInput() {
    renameFileToInput('/rename-but1-copy-to-input');
}

function renameBut2ToInput() {
    renameFileToInput('/rename-but2-copy-to-input');
}

function renameBut3ToInput() {
    renameFileToInput('/rename-but3-copy-to-input');
}

function renameBut4ToInput() {
    renameFileToInput('/rename-but4-copy-to-input');
}

function renameBut5ToInput() {
    renameFileToInput('/rename-but5-copy-to-input');
}

  // Fonction pour envoyer une requête POST au serveur
  function envoyerMot(mot) {
    fetch('http://localhost:3000/write', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ word: mot })
    })
    .then(response => response.text())
    .then(data => {
        console.log(data);
    })
    .catch(error => {
        console.error('Erreur lors de l\'envoi de la requête:', error);
    });
}

// Fonction pour envoyer une requête POST au serveur
function envoyerMot1(mot) {
    fetch('http://localhost:3000/write1', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ word: mot })
    })
    .then(response => response.text())
    .then(data => {
        console.log(data);
    })
    .catch(error => {
        console.error('Erreur lors de l\'envoi de la requête:', error);
    });
}

// Fonction pour envoyer une requête POST au serveur
function envoyerMot2(mot) {
    fetch('http://localhost:3000/write2', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({ word: mot })
    })
    .then(response => response.text())
    .then(data => {
        console.log(data);
    })
    .catch(error => {
        console.error('Erreur lors de l\'envoi de la requête:', error);
    });
}





   


