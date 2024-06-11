const express = require('express');
const pl = require('tau-prolog');
const fs = require('fs').promises;
const path = require('path');
const browserSync = require('browser-sync');
const { exec } = require('child_process');

const app = express();
const port = process.env.PORT || 3000;
const session = pl.create(100000000);


// Middleware pour servir les fichiers statiques depuis le répertoire 'public'
app.use(express.static(path.join(__dirname, 'public')));
app.use(express.json());
app.use(express.urlencoded({ extended: true }));



// Configuration de BrowserSync
const bs = browserSync.create();
bs.init({
    proxy: "http://localhost:3000",
    files: ["public/**/*.*"],
    open: false,
    notify: false
});

// Charger le programme Prolog depuis le fichier code4.pl
const loadPrologProgram = async () => {
    try {
        const filePath = path.join(__dirname, 'code4.pl');
        const data = await fs.readFile(filePath, 'utf8');
        return data;
    } catch (err) {
        throw new Error('Erreur lors de la lecture du fichier Prolog : ' + err.message);
    }
};


// Charger le programme Prolog au démarrage du serveur
app.listen(port, async () => {
    try {
        const program = await loadPrologProgram();
        await session.consult(program);
        console.log('Programme Prolog chargé avec succès');
        console.log('Serveur démarré sur le port', port);
        
    } catch (error) {
        console.error('Erreur lors du chargement du programme Prolog :', error);
        process.exit(1); // Arrêter le serveur en cas d'erreur critique
    }
});

const executePrologQuery = async (queryText) => {
    return new Promise((resolve, reject) => {
        session.query(queryText, {
            success: () => {
                session.answers(answer => {
                    const formattedAnswer = session.format_answer(answer);
                    resolve(formattedAnswer);
                });
            },
            error: (err) => {
                reject(new Error('Erreur lors de l\'exécution de la requête Prolog : ' + err.message));
            }
        });
    });
};


// Route POST pour exécuter les requêtes Prolog
app.post('/execute-prolog', async (req, res) => {
    const queryText = req.body.query;
    console.log('Requête Prolog reçue :', queryText);

    try {
        // Exécuter la requête Prolog
        const result = await executePrologQuery(queryText);
        console.log('Résultat de la requête Prolog :', result);
        
         // Écrire le résultat dans le fichier fresult.txt
         await fs.writeFile('fresult.txt', result);
        // Envoyer les résultats au format JSON dans la réponse
        res.json({ results: [result] }); // Envelopper le résultat dans un tableau
    } catch (error) {
        console.error('Erreur lors de l\'exécution de la requête Prolog :', error);
        res.status(500).json({ error: 'Erreur lors de l\'exécution de la requête Prolog' });
    }
});



// Route pour copier but1.txt en but1_copy.txt et renommer but1_copy.txt en input.txt
app.post('/rename-but1-copy-to-input', async (req, res) => {
    try {
        
        const inputFile = 'input.txt';
        const but1File = 'but1.txt';
        const but1CopyFile = 'but1_copy.txt';
        const outputFile = 'output.txt';

        // Copier but1.txt en but1_copy.txt
        await fs.copyFile(but1File, but1CopyFile);

        // Supprimer input.txt s'il existe
        await fs.unlink(inputFile).catch(err => {
            if (err.code !== 'ENOENT') {
                throw err;
            }
        });

        // Renommer but1Copy.txt en input.txt
        await fs.rename(but1CopyFile, inputFile);

         // Supprimer le contenu de output.txt
         await fs.writeFile(outputFile, ''); // Écrire une chaîne vide pour vider le fichier

        res.send('Renommage réussi et output.txt vidé');
    } catch (error) {
        console.error('Erreur lors du renommage du fichier:', error);
        res.status(500).send('Erreur lors du renommage du fichier');
    }
});


// Route pour copier but2.txt en but2_copy.txt et renommer but2_copy.txt en input.txt
app.post('/rename-but2-copy-to-input', async (req, res) => {
    try {
        const inputFile = 'input.txt';
        const but2File = 'but2.txt';
        const but2CopyFile = 'but2_copy.txt';
        const outputFile = 'output.txt';

        // Copier but2.txt en but2_copy.txt
        await fs.copyFile(but2File, but2CopyFile);

        // Supprimer input.txt s'il existe
        await fs.unlink(inputFile).catch(err => {
            if (err.code !== 'ENOENT') {
                throw err;
            }
        });

        // Renommer but2_copy.txt en input.txt
        await fs.rename(but2CopyFile, inputFile);

        // Supprimer le contenu de output.txt
        await fs.writeFile(outputFile, ''); // Écrire une chaîne vide pour vider le fichier

        res.send('Renommage réussi et output.txt vidé');
    } catch (error) {
        console.error('Erreur lors du renommage du fichier:', error);
        res.status(500).send('Erreur lors du renommage du fichier');
    }
});




// Route pour copier but3.txt en but3_copy.txt et renommer but3_copy.txt en input.txt
app.post('/rename-but3-copy-to-input', async (req, res) => {
    try {
        const inputFile = 'input.txt';
        const but3File = 'but3.txt';
        const but3CopyFile = 'but3_copy.txt';
        const outputFile = 'output.txt';

        // Copier but3.txt en but3_copy.txt
        await fs.copyFile(but3File, but3CopyFile);

        // Supprimer input.txt s'il existe
        await fs.unlink(inputFile).catch(err => {
            if (err.code !== 'ENOENT') {
                throw err;
            }
        });

        // Renommer but3_copy.txt en input.txt
        await fs.rename(but3CopyFile, inputFile);

        // Supprimer le contenu de output.txt
        await fs.writeFile(outputFile, ''); // Écrire une chaîne vide pour vider le fichier

        res.send('Renommage réussi et output.txt vidé');
    } catch (error) {
        console.error('Erreur lors du renommage du fichier:', error);
        res.status(500).send('Erreur lors du renommage du fichier');
    }
});


// Route pour copier but4.txt en but4_copy.txt et renommer but4_copy.txt en input.txt
app.post('/rename-but4-copy-to-input', async (req, res) => {
    try {
        const inputFile = 'input.txt';
        const but4File = 'but4.txt';
        const but4CopyFile = 'but4_copy.txt';
        const outputFile = 'output.txt';

        // Copier but4.txt en but4_copy.txt
        await fs.copyFile(but4File, but4CopyFile);

        // Supprimer input.txt s'il existe
        await fs.unlink(inputFile).catch(err => {
            if (err.code !== 'ENOENT') {
                throw err;
            }
        });

        // Renommer but4_copy.txt en input.txt
        await fs.rename(but4CopyFile, inputFile);

        // Supprimer le contenu de output.txt
        await fs.writeFile(outputFile, ''); // Écrire une chaîne vide pour vider le fichier

        res.send('Renommage réussi et output.txt vidé');
    } catch (error) {
        console.error('Erreur lors du renommage du fichier:', error);
        res.status(500).send('Erreur lors du renommage du fichier');
    }
});





// Route pour copier but5.txt en but5_copy.txt et renommer but5_copy.txt en input.txt
app.post('/rename-but5-copy-to-input', async (req, res) => {
    try {
        const inputFile = 'input.txt';
        const but5File = 'but5.txt';
        const but5CopyFile = 'but5_copy.txt';
        const outputFile = 'output.txt';

        // Copier but5.txt en but5_copy.txt
        await fs.copyFile(but5File, but5CopyFile);

        // Supprimer input.txt s'il existe
        await fs.unlink(inputFile).catch(err => {
            if (err.code !== 'ENOENT') {
                throw err;
            }
        });

        // Renommer but5_copy.txt en input.txt
        await fs.rename(but5CopyFile, inputFile);

        // Supprimer le contenu de output.txt
        await fs.writeFile(outputFile, ''); // Écrire une chaîne vide pour vider le fichier

        res.send('Renommage réussi et output.txt vidé');
    } catch (error) {
        console.error('Erreur lors du renommage du fichier:', error);
        res.status(500).send('Erreur lors du renommage du fichier');
    }
});


// route GET pour récupérer le contenu du fichier output.txt
app.get('/get-output', async (req, res) => {
    try {
        const output = await fs.readFile('output.txt', 'utf-8');
        res.send(output);
    } catch (error) {
        console.error('Erreur lors de la lecture du fichier output.txt :', error);
        res.status(500).json({ error: 'Erreur lors de la lecture du fichier output.txt' });
    }
});







// Route pour servir le contenu de code_exemples.pl
app.get('/get-code-exemples', async (req, res) => {
    try {
        const codeExemplesPath = path.join(__dirname, 'code_exemples.pl');
        const codeExemplesContent = await fs.readFile(codeExemplesPath, 'utf8');
        res.send(codeExemplesContent);
    } catch (error) {
        console.error('Erreur lors de la lecture de code_exemples.pl :', error);
        res.status(500).json({ error: 'Erreur lors de la lecture de code_exemples.pl' });
    }
});

// Route GET pour afficher le script combiné à partir de input.txt et getg.txt
app.get('/afficher-script', async (req, res) => {
    try {
        // Lire le contenu de input.txt
        const inputPath = path.join(__dirname, 'input.txt');
        const inputContent = await fs.readFile(inputPath, 'utf8');

        // Lire le contenu de getg.txt
        const getgPath = path.join(__dirname, 'getg.txt');
        const getgContent = await fs.readFile(getgPath, 'utf8');

        // Combiner les deux contenus
        const scriptContent = inputContent.split('\n');
        scriptContent.splice(3, 0, getgContent);
        const combinedScript = scriptContent.join('\n');

        // Envoyer le contenu combiné en réponse à la demande
        res.send(combinedScript);
    } catch (error) {
        console.error('Erreur lors de la lecture des fichiers ou de la combinaison des contenus :', error);
        res.status(500).send('Erreur lors de la récupération du script combiné');
    }
});















// Route pour écrire les mots dans le fichier getg.txt
app.post('/write', async (req, res) => {
    const { word } = req.body;

    try {
        // Effacer le contenu de getg.txt avant d'écrire les nouvelles données
       await fs.writeFile('getg.txt', '', { flag: 'w' });

        // Écrire le mot dans getg.txt
        await fs.appendFile('getg.txt', word + '\n');

        res.send('Écriture réussie');
    } catch (err) {
        console.error('Erreur lors de l\'écriture dans le fichier:', err);
        res.status(500).send('Erreur lors de l\'écriture dans le fichier');
    }
});




// Route pour trouver, effacer et écrire les mots dans le fichier input.txt
app.post('/write1', async (req, res) => {
    const { word } = req.body;

    try {
        // Vérifier si le fichier input.txt existe
        const fileExists = await fs.access('input.txt').then(() => true).catch(() => false);
        
        if (fileExists) {
            // Effacer le contenu de input.txt avant d'écrire les nouvelles données
            await fs.writeFile('input.txt', '', { flag: 'w' });
            await fs.writeFile('output.txt', '', { flag: 'w' });
        }

        // Écrire le mot dans input.txt
        await fs.appendFile('input.txt', word + '\n');

        res.send('Écriture réussie');
    } catch (err) {
        console.error('Erreur lors de l\'écriture dans le fichier:', err);
        res.status(500).send('Erreur lors de l\'écriture dans le fichier');
    }
});



// Route pour seulement écrire les mots dans le fichier input.txt
app.post('/write2', async (req, res) => {
    const { word } = req.body;

    try {
        // Écrire le mot dans input.txt
        await fs.appendFile('input.txt', word + '\n');

        res.send('Écriture réussie');
    } catch (err) {
        console.error('Erreur lors de l\'écriture dans le fichier:', err);
        res.status(500).send('Erreur lors de l\'écriture dans le fichier');
    }
});



// Définir la route pour tree.html
app.get('/tree.html', (req, res) => {
    res.sendFile(path.join(__dirname, 'public', 'tree.html'));
  });



// Route POST pour effacer le contenu du fichier output.txt
app.post('/effacer-output', async (req, res) => {
    try {
        // Effacer le contenu du fichier output.txt
        const outputFile = 'output.txt';
        await fs.writeFile(outputFile, '');
        res.sendStatus(200); // Réponse indiquant que l'effacement a réussi
    } catch (error) {
        console.error('Erreur lors de l\'effacement du fichier output.txt :', error);
        res.status(500).json({ error: 'Erreur lors de l\'effacement du fichier output.txt' });
    }
});






