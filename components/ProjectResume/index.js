import React from "react";

const ProjectResume = ({ dates, type, location, position, bullets }) => {
  const [bulletsLocal, setBulletsLocal] = React.useState(bullets.split(","));

  return (
    <div className="mt-5 w-full flex mob:flex-col desktop:flex-row justify-between">
      <div className="text-lg w-2/5">
        <h2>{dates}</h2>
        <h3 className="text-sm opacity-50">{location}</h3>
        <h3 className="text-sm opacity-50">{type}</h3>
      </div>
      <div className="w-3/5">
        <h2 className="text-lg font-bold">{position}</h2>
        {bulletsLocal && bulletsLocal.length > 0 && (
          <ul className="list-disc">
            {bulletsLocal
              .join("") // Convert the array to a single string if needed
              .split("#") // Split the string by the "#" character
              .filter((bullet) => bullet.trim().length > 0) // Remove empty items
              .map((bullet, index) => (
                <li key={index} className="text-sm my-1 opacity-70">
                  {bullet.trim()} {/* Remove any extra spaces */}
                </li>
              ))}
          </ul>
)}

      </div>
    </div>
  );
};

export default ProjectResume;
